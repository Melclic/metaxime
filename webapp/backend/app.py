import enum
import json
import os
import queue
import subprocess
import threading
import uuid
from contextlib import asynccontextmanager
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

from metaxime.utils import convert_depiction

# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------

SCRIPT_DIR: Path = Path(__file__).resolve().parent

WORKFLOW_SCRIPT: Path = (SCRIPT_DIR / "../nextflow/complete_workflow.nf").resolve()
WORKFLOW_CONFIG: Path = (SCRIPT_DIR / "../nextflow/nextflow.config").resolve()

RUN_DB_PATH: Path = SCRIPT_DIR / "run_db.json"
DEFAULT_RUN_ROOT: Path = SCRIPT_DIR / "run_files"

# ---------------------------------------------------------------------------
# Job models
# ---------------------------------------------------------------------------


class JobStatus(str, enum.Enum):
    """Enumeration of possible job states."""

    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


class JobCreateRequest(BaseModel):
    """Request body for creation of a new workflow job.

    A per-job folder will be created under base_output_dir:
        <base_output_dir>/<job_id>/
          _work_dir/   -> Nextflow -work-dir
          _log/        -> Nextflow log file
          _output/     -> workflow outputs (used as --output_folder)

    Attributes:
        model_file:
            Path to the model file passed to --model_file.
        target_inchi:
            InChI string passed to --source_inchi.
        base_output_dir:
            Root directory where per job folders will be created.
            If omitted, a default "run_files" folder next to app.py is used.
        work_dir:
            Optional working directory used as subprocess cwd.
        rules_file:
            Optional path for --rules_file.
        std_mode:
            Optional string for --std_mode.
        max_steps:
            Optional integer for --max_steps.
        topx:
            Optional integer for --topx.
        accept_partial_results:
            Boolean flag for --accept_partial_results.
        diameters:
            Optional comma separated diameters string for --diameters.
        rule_type:
            Optional rule type for --rule_type.
        ram_limit:
            Optional integer ram limit in GB for --ram_limit.
        source_comp:
            Optional source compartment id for --source_comp.
        target_comp:
            Optional target compartment id for --target_comp.
        use_inchikey2:
            Boolean flag for --use_inchikey2.
        find_all_parentless:
            Boolean flag for --find_all_parentless.
    """

    model_file: str
    target_inchi: str

    base_output_dir: Optional[str] = None
    work_dir: Optional[str] = None

    rules_file: Optional[str] = None
    std_mode: Optional[str] = None
    max_steps: Optional[int] = None
    topx: Optional[int] = None
    accept_partial_results: Optional[bool] = None
    diameters: Optional[str] = None
    rule_type: Optional[str] = None
    ram_limit: Optional[int] = None

    source_comp: Optional[str] = None
    target_comp: Optional[str] = None
    use_inchikey2: Optional[bool] = None
    find_all_parentless: Optional[bool] = None


class JobInfo(BaseModel):
    """Representation of a workflow job and its state.

    Attributes:
        id:
            Unique identifier of the job.
        status:
            Current status of the job.
        created_at:
            Time when the job was created.
        started_at:
            Time when the job was picked up by the worker.
        finished_at:
            Time when the job finished.
        work_dir:
            Working directory used as subprocess cwd.
        job_dir:
            Root directory created for this job.
        nxf_work_dir:
            Directory used for Nextflow -work-dir.
        log_dir:
            Directory for Nextflow log files.
        log_file:
            Path to the Nextflow log file (-log).
        output_folder:
            Directory used for --output_folder.
        payload:
            Original request parameters plus derived values.
        command_preview:
            Shell quoted preview of the executed command.
        exit_code:
            Process exit code if the job has finished.
        error_message:
            Short error message when the job failed.
        stdout:
            Captured standard output of the workflow.
        stderr:
            Captured standard error of the workflow.
    """

    id: str
    status: JobStatus
    created_at: datetime
    started_at: Optional[datetime] = None
    finished_at: Optional[datetime] = None

    work_dir: Optional[str] = None
    job_dir: Optional[str] = None
    nxf_work_dir: Optional[str] = None
    log_dir: Optional[str] = None
    log_file: Optional[str] = None
    output_folder: Optional[str] = None

    payload: Dict[str, Any]
    command_preview: Optional[str] = None

    exit_code: Optional[int] = None
    error_message: Optional[str] = None
    stdout: Optional[str] = None
    stderr: Optional[str] = None


class JobSummary(BaseModel):
    """Light version of a job record for list endpoints."""

    id: str
    status: JobStatus
    created_at: datetime
    started_at: Optional[datetime] = None
    finished_at: Optional[datetime] = None
    target_inchi: str = None
    target_smiles: str = None


# ---------------------------------------------------------------------------
# In memory store and synchronization with JSON DB
# ---------------------------------------------------------------------------

job_store: Dict[str, JobInfo] = {}
job_store_lock = threading.Lock()

job_queue: "queue.Queue[str]" = queue.Queue()
worker_thread_started: bool = False
worker_thread_lock = threading.Lock()

run_db_lock = threading.Lock()


def shell_join(args: List[str]) -> str:
    """Return a shell escaped representation of an argument list.

    Args:
        args:
            Command argument list.

    Returns:
        str:
            Shell quoted string useful for logging or debugging.
    """
    import shlex

    return " ".join(shlex.quote(a) for a in args)


def load_job_store_from_db() -> None:
    """Load persisted jobs from RUN_DB_PATH into job_store.

    The JSON file is expected to contain a list of JobInfo dictionaries.
    Any malformed entries are skipped. Jobs in RUNNING or QUEUED state
    at load time are marked as FAILED since they cannot still be running
    after a restart.
    """
    global job_store

    if not RUN_DB_PATH.exists():
        return

    with run_db_lock:
        try:
            content = RUN_DB_PATH.read_text(encoding="utf8").strip()
            if not content:
                return
            data = json.loads(content)
        except Exception:
            return

    if not isinstance(data, list):
        return

    loaded: Dict[str, JobInfo] = {}
    for entry in data:
        if not isinstance(entry, dict):
            continue
        try:
            job = JobInfo(**entry)
        except Exception:
            continue

        if job.status in {JobStatus.RUNNING, JobStatus.QUEUED}:
            job.status = JobStatus.FAILED
            if job.error_message is None:
                job.error_message = "Server restarted while job was in progress"
            if job.finished_at is None:
                job.finished_at = datetime.utcnow()

        loaded[job.id] = job

    with job_store_lock:
        job_store = loaded


def persist_job_store() -> None:
    """Persist the entire job_store to RUN_DB_PATH atomically.

    This function should be called after any mutation of job_store to keep
    the JSON database and in memory store in sync.
    """
    with job_store_lock:
        data = [job.model_dump() for job in job_store.values()]

    serialized = json.dumps(data, default=str, indent=2)

    with run_db_lock:
        tmp_path = RUN_DB_PATH.with_suffix(".tmp")
        tmp_path.write_text(serialized, encoding="utf8")
        tmp_path.replace(RUN_DB_PATH)


# ---------------------------------------------------------------------------
# Command construction and worker logic
# ---------------------------------------------------------------------------


def build_nextflow_command(job: JobInfo) -> List[str]:
    """Build the Nextflow command argument list for a given job.

    Args:
        job:
            JobInfo instance that contains the workflow payload.

    Returns:
        List[str]:
            Argument vector suitable for subprocess execution.
    """
    p: Dict[str, Any] = job.payload

    model_file = str(p["model_file"])
    target_inchi = str(p["target_inchi"])
    output_folder = str(p["output_folder"])
    nxf_work_dir = str(p["nxf_work_dir"])
    log_file = str(p["log_file"])

    args: List[str] = [
        "nextflow",
        "run",
        str(WORKFLOW_SCRIPT),
        "--model_file",
        model_file,
        "--source_inchi",
        target_inchi,
        "--output_folder",
        output_folder,
        "-work-dir",
        nxf_work_dir,
        "-c",
        str(WORKFLOW_CONFIG),
        "-ansi-log",
        "false",
        "-log",
        log_file,
    ]

    if p.get("rules_file") is not None:
        args.extend(["--rules_file", str(p["rules_file"])])

    if p.get("std_mode") is not None:
        args.extend(["--std_mode", str(p["std_mode"])])

    if p.get("max_steps") is not None:
        args.extend(["--max_steps", str(p["max_steps"])])

    if p.get("topx") is not None:
        args.extend(["--topx", str(p["topx"])])

    if p.get("accept_partial_results"):
        args.append("--accept_partial_results")

    if p.get("diameters") is not None:
        args.extend(["--diameters", str(p["diameters"])])

    if p.get("rule_type") is not None:
        args.extend(["--rule_type", str(p["rule_type"])])

    if p.get("ram_limit") is not None:
        args.extend(["--ram_limit", str(p["ram_limit"])])

    if p.get("source_comp") is not None:
        args.extend(["--source_comp", str(p["source_comp"])])

    if p.get("target_comp") is not None:
        args.extend(["--target_comp", str(p["target_comp"])])

    if p.get("use_inchikey2"):
        args.append("--use_inchikey2")

    if p.get("find_all_parentless"):
        args.append("--find_all_parentless")

    return args


def worker_loop() -> None:
    """Background worker that runs queued jobs one at a time.

    This function is intended to run in a dedicated daemon thread.
    It waits for job identifiers from the queue, updates job state,
    executes Nextflow, captures output, and persists state changes
    to the JSON database.
    """
    while True:
        job_id = job_queue.get()
        try:
            with job_store_lock:
                job = job_store.get(job_id)

            if job is None:
                continue

            with job_store_lock:
                job.status = JobStatus.RUNNING
                job.started_at = datetime.utcnow()
                job_store[job_id] = job
            persist_job_store()

            args = build_nextflow_command(job)
            args = [str(a) for a in args]

            with job_store_lock:
                job.command_preview = shell_join(args)
                job_store[job_id] = job
            persist_job_store()

            try:
                process = subprocess.Popen(
                    args,
                    cwd=job.job_dir,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                stdout, stderr = process.communicate()
                exit_code = process.returncode

                with job_store_lock:
                    job.stdout = stdout
                    job.stderr = stderr
                    job.exit_code = exit_code
                    job.finished_at = datetime.utcnow()
                    if exit_code == 0:
                        job.status = JobStatus.COMPLETED
                    else:
                        job.status = JobStatus.FAILED
                        job.error_message = f"Nextflow exited with code {exit_code}"
                    job_store[job_id] = job
                persist_job_store()
            except Exception as exc:
                with job_store_lock:
                    job.status = JobStatus.FAILED
                    job.error_message = str(exc)
                    job.finished_at = datetime.utcnow()
                    job_store[job_id] = job
                persist_job_store()
        finally:
            job_queue.task_done()


def ensure_worker_thread_started() -> None:
    """Start the worker thread once for the application lifetime."""
    global worker_thread_started
    with worker_thread_lock:
        if not worker_thread_started:
            thread = threading.Thread(target=worker_loop, daemon=True)
            thread.start()
            worker_thread_started = True


# ---------------------------------------------------------------------------
# FastAPI app and endpoints
# ---------------------------------------------------------------------------


@asynccontextmanager
async def lifespan(app: FastAPI):
    """FastAPI lifespan context.

    Loads the job store from disk, starts the worker, and persists
    state once more on shutdown.
    """
    DEFAULT_RUN_ROOT.mkdir(parents=True, exist_ok=True)
    load_job_store_from_db()
    ensure_worker_thread_started()
    yield
    persist_job_store()


app = FastAPI(
    title="Nextflow Workflow Queue API",
    version="1.0.0",
    lifespan=lifespan,
)

from fastapi.middleware.cors import CORSMiddleware

# Allow the React dev server origin
origins = [
    "http://localhost:5173",
    "http://127.0.0.1:5173",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,          # or ["*"] for everything in dev
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.post("/jobs", response_model=JobInfo)
def create_job(request: JobCreateRequest) -> JobInfo:
    """Create a new workflow job and enqueue it.

    A per job directory structure is created and added to the payload.
    The job is stored in memory and persisted to the JSON database.

    Args:
        request:
            Workflow request payload.

    Returns:
        JobInfo:
            Full job description including its initial state.
    """
    payload: Dict[str, Any] = request.model_dump()
    run_cwd: Optional[str] = payload.pop("work_dir", None)
    base_output_dir_str: Optional[str] = payload.pop("base_output_dir", None)

    base_output_dir = (
        Path(base_output_dir_str).resolve()
        if base_output_dir_str
        else DEFAULT_RUN_ROOT
    )

    job_id: str = str(uuid.uuid4())

    job_dir_path = base_output_dir / job_id
    nxf_work_dir_path = job_dir_path / "_work_dir"
    log_dir_path = job_dir_path / "_log"
    output_folder_path = job_dir_path / "_output"
    log_file_path = log_dir_path / "nextflow.log"

    nxf_work_dir_path.mkdir(parents=True, exist_ok=True)
    log_dir_path.mkdir(parents=True, exist_ok=True)
    output_folder_path.mkdir(parents=True, exist_ok=True)

    job_dir = str(job_dir_path)
    nxf_work_dir = str(nxf_work_dir_path)
    log_dir = str(log_dir_path)
    log_file = str(log_file_path)
    output_folder = str(output_folder_path)

    payload["model_file"] = str(payload["model_file"])
    payload["target_inchi"] = str(payload["target_inchi"])
    payload["target_smiles"] = convert_depiction(str(payload["target_inchi"]), 'inchi', 'smiles')
    payload["job_dir"] = job_dir
    payload["nxf_work_dir"] = nxf_work_dir
    payload["log_dir"] = log_dir
    payload["log_file"] = log_file
    payload["output_folder"] = output_folder

    job = JobInfo(
        id=job_id,
        status=JobStatus.QUEUED,
        created_at=datetime.utcnow(),
        work_dir=run_cwd,
        job_dir=job_dir,
        nxf_work_dir=nxf_work_dir,
        log_dir=log_dir,
        log_file=log_file,
        output_folder=output_folder,
        payload=payload,
    )

    with job_store_lock:
        job_store[job_id] = job
    persist_job_store()

    job_queue.put(job_id)

    return job


@app.get("/jobs", response_model=List[JobSummary])
def list_jobs() -> List[JobSummary]:
    """List all jobs with basic information.

    Returns:
        List[JobSummary]:
            Sorted by creation time (newest first).
    """
    with job_store_lock:
        jobs = sorted(
            job_store.values(),
            key=lambda j: j.created_at,
            reverse=True,
        )

    return [
        JobSummary(
            id=j.id,
            status=j.status,
            created_at=j.created_at,
            started_at=j.started_at,
            finished_at=j.finished_at,
            target_inchi=j.payload['target_inchi'],
            target_smiles=j.payload['target_smiles'],
        )
        for j in jobs
    ]


@app.get("/jobs/running", response_model=List[JobSummary])
def list_running_jobs() -> List[JobSummary]:
    """List jobs that are currently running."""
    with job_store_lock:
        jobs = [
            j
            for j in job_store.values()
            if j.status == JobStatus.RUNNING
        ]

    jobs = sorted(jobs, key=lambda j: j.created_at, reverse=True)

    return [
        JobSummary(
            id=j.id,
            status=j.status,
            created_at=j.created_at,
            started_at=j.started_at,
            finished_at=j.finished_at,
            target_inchi=j.payload['target_inchi'],
            target_smiles=j.payload['target_smiles'],
        )
        for j in jobs
    ]


@app.get("/jobs/completed", response_model=List[JobSummary])
def list_completed_jobs() -> List[JobSummary]:
    """List jobs that have completed successfully."""
    with job_store_lock:
        jobs = [
            j
            for j in job_store.values()
            if j.status == JobStatus.COMPLETED
        ]

    jobs = sorted(jobs, key=lambda j: j.created_at, reverse=True)

    return [
        JobSummary(
            id=j.id,
            status=j.status,
            created_at=j.created_at,
            started_at=j.started_at,
            finished_at=j.finished_at,
            target_inchi=j.payload['target_inchi'],
            target_smiles=j.payload['target_smiles'],
        )
        for j in jobs
    ]


@app.get("/jobs/failed", response_model=List[JobSummary])
def list_failed_jobs() -> List[JobSummary]:
    """List jobs that have failed."""
    with job_store_lock:
        jobs = [
            j
            for j in job_store.values()
            if j.status == JobStatus.FAILED
        ]

    jobs = sorted(jobs, key=lambda j: j.created_at, reverse=True)

    return [
        JobSummary(
            id=j.id,
            status=j.status,
            created_at=j.created_at,
            started_at=j.started_at,
            finished_at=j.finished_at,
            target_inchi=j.payload['target_inchi'],
            target_smiles=j.payload['target_smiles'],
        )
        for j in jobs
    ]


@app.get("/jobs/{job_id}", response_model=JobInfo)
def get_job(job_id: str) -> JobInfo:
    """Retrieve full information for a single job.

    Args:
        job_id:
            Identifier of the job.

    Returns:
        JobInfo:
            Full job record.

    Raises:
        HTTPException:
            If the job does not exist.
    """
    with job_store_lock:
        job = job_store.get(job_id)

    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")

    return job


@app.delete("/jobs/{job_id}", response_model=JobSummary)
def delete_job(job_id: str) -> JobSummary:
    """Remove a finished job from the in memory store.

    Only COMPLETED or FAILED jobs can be deleted. The
    corresponding entry is removed from the JSON database
    via persist_job_store.

    Args:
        job_id:
            Identifier of the job to delete.

    Returns:
        JobSummary:
            Summary of the removed job.

    Raises:
        HTTPException:
            If the job does not exist or is not yet finished.
    """
    with job_store_lock:
        job = job_store.get(job_id)
        if job is None:
            raise HTTPException(status_code=404, detail="Job not found")
        if job.status in {JobStatus.QUEUED, JobStatus.RUNNING}:
            raise HTTPException(
                status_code=400,
                detail="Cannot delete a queued or running job",
            )
        summary = JobSummary(
            id=job.id,
            status=job.status,
            created_at=job.created_at,
            started_at=job.started_at,
            finished_at=job.finished_at,
            target_inchi=j.payload['target_inchi'],
            target_smiles=j.payload['target_smiles'],
        )
        del job_store[job_id]

    persist_job_store()
    return summary


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(
        app=app,
        host="0.0.0.0",
        port=8000,
        reload=False,
    )