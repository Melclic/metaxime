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
from typing import Any
from typing import Dict
from typing import List
from typing import Optional

from fastapi import FastAPI
from fastapi import HTTPException
from pydantic import BaseModel

SCRIPT_DIR = Path(__file__).resolve().parent

RUN_DB_PATH: Path = SCRIPT_DIR / "run_db.json"
RUN_FILES_PATH: Path = SCRIPT_DIR / "run_files/"
RUN_FILES_PATH.mkdir(parents=True, exist_ok=True)

WORKFLOW_SCRIPT = SCRIPT_DIR / "../nextflow/complete_workflow.nf"
WORKFLOW_CONFIG = SCRIPT_DIR / "../nextflow/nextflow.config"

job_store: Dict[str, "JobInfo"] = {}
job_store_lock = threading.Lock()

job_queue: "queue.Queue[str]" = queue.Queue()
worker_thread_started: bool = False
worker_thread_lock = threading.Lock()

run_db_lock = threading.Lock()


class JobStatus(str, enum.Enum):
    """Enumeration of possible job states."""

    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


class JobCreateRequest(BaseModel):
    """Request body for creation of a new Nextflow workflow job.

    The required parameters map to the base argument list.
    Optional fields correspond to extra workflow flags and
    are only passed when provided.

    Attributes:
        model_file:
            Path to the model file passed to --model_file.
        target_inchi:
            InChI string passed to --source_inchi.
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
    """Representation of a Nextflow workflow job and its state.

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
        job_dir:
            Root directory created for this job.
        nxf_work_dir:
            Directory used for Nextflow -work-dir.
        log_dir:
            Directory where Nextflow log is written.
        output_tar_path:
            Path used for --out_tar.
        payload:
            Original request parameters for the workflow plus derived values.
        command_preview:
            Optional string preview of the command arguments.
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

    job_dir: Optional[Path] = None
    nxf_work_dir: Optional[Path] = None
    output_folder: Optional[Path] = None

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


def build_nextflow_command(
    job: JobInfo,
) -> List[str]:
    """Build the Nextflow argument list for a given job."""

    p: Dict[str, Any] = job.payload

    model_file: str = p["model_file"]
    target_inchi: str = p["target_inchi"]
    nxf_work_dir: str = p["nxf_work_dir"]
    output_folder: str = p["output_folder"]

    # if not target_inchi[0]=='"' and not target_inchi[-1]=='"':
    #     target_inchi_quoted = f'"{target_inchi}"'
    # else:
    #     target_inchi_quoted = target_inchi

    args: List[str] = [
        "nextflow",
        "run",
        str(WORKFLOW_SCRIPT.resolve()),
        "--model_file",
        model_file,
        "--source_inchi",
        target_inchi,
        "--output_folder",
        str(output_folder.resolve()),
        "-work-dir",
        str(nxf_work_dir.resolve()),
        "-c",
        str(WORKFLOW_CONFIG.resolve()),
        "-ansi-log",
        "false",
    ]

    if p.get("rules_file") is not None:
        args.extend(
            [
                "--rules_file",
                str(p["rules_file"]),
            ],
        )

    if p.get("std_mode") is not None:
        args.extend(
            [
                "--std_mode",
                str(p["std_mode"]),
            ],
        )

    if p.get("max_steps") is not None:
        args.extend(
            [
                "--max_steps",
                str(p["max_steps"]),
            ],
        )

    if p.get("topx") is not None:
        args.extend(
            [
                "--topx",
                str(p["topx"]),
            ],
        )

    if p.get("accept_partial_results"):
        args.append("--accept_partial_results")

    if p.get("diameters") is not None:
        args.extend(
            [
                "--diameters",
                str(p["diameters"]),
            ],
        )

    if p.get("rule_type") is not None:
        args.extend(
            [
                "--rule_type",
                str(p["rule_type"]),
            ],
        )

    if p.get("ram_limit") is not None:
        args.extend(
            [
                "--ram_limit",
                str(p["ram_limit"]),
            ],
        )

    if p.get("source_comp") is not None:
        args.extend(
            [
                "--source_comp",
                str(p["source_comp"]),
            ],
        )

    if p.get("target_comp") is not None:
        args.extend(
            [
                "--target_comp",
                str(p["target_comp"]),
            ],
        )

    if p.get("use_inchikey2"):
        args.append("--use_inchikey2")

    if p.get("find_all_parentless"):
        args.append("--find_all_parentless")

    return args


def append_run_to_db(
    job: JobInfo,
) -> None:
    """Append a finished job record to the JSON run database."""

    with run_db_lock:
        data: List[Dict[str, Any]] = []

        if RUN_DB_PATH.exists():
            try:
                with RUN_DB_PATH.open(
                    mode="r",
                    encoding="utf8",
                ) as handle:
                    content = handle.read().strip()
                    if content:
                        data = json.loads(content)
            except json.JSONDecodeError:
                data = []

        data.append(job.model_dump())

        with RUN_DB_PATH.open(
            mode="w",
            encoding="utf8",
        ) as handle:
            json.dump(
                obj=data,
                fp=handle,
                default=str,
                indent=2,
            )


def worker_loop() -> None:
    """Background worker that runs queued jobs one at a time."""

    while True:
        job_id = job_queue.get()
        try:
            with job_store_lock:
                job = job_store.get(job_id)

            if job is None:
                continue

            with job_store_lock:
                job.status = JobStatus.RUNNING
                job.started_at = datetime.now()
                job_store[job_id] = job

            args: List[str] = build_nextflow_command(job)
            with job_store_lock:
                job.command_preview = " ".join(args)
                job_store[job_id] = job

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
                    job.finished_at = datetime.now()
                    if exit_code == 0:
                        job.status = JobStatus.COMPLETED
                    else:
                        job.status = JobStatus.FAILED
                        job.error_message = (
                            f"Nextflow exited with code {exit_code}"
                        )
                    job_store[job_id] = job

                if job.status in {
                    JobStatus.COMPLETED,
                    JobStatus.FAILED,
                }:
                    try:
                        append_run_to_db(job)
                    except Exception:
                        pass
            except Exception as exc:
                with job_store_lock:
                    job.status = JobStatus.FAILED
                    job.error_message = str(exc)
                    job.finished_at = datetime.now()
                    job_store[job_id] = job

                try:
                    append_run_to_db(job)
                except Exception:
                    pass
        finally:
            job_queue.task_done()


def ensure_worker_thread_started() -> None:
    """Start the worker thread once for the application lifetime."""

    global worker_thread_started
    with worker_thread_lock:
        if not worker_thread_started:
            thread = threading.Thread(
                target=worker_loop,
                daemon=True,
            )
            thread.start()
            worker_thread_started = True


@asynccontextmanager
async def lifespan(
    app: FastAPI,
):
    """FastAPI lifespan context."""

    ensure_worker_thread_started()
    yield


app = FastAPI(
    title="Nextflow Workflow Queue API",
    version="1.0.0",
    lifespan=lifespan,
)


@app.post(
    "/jobs",
    response_model=JobInfo,
)
def create_job(
    request: JobCreateRequest,
) -> JobInfo:
    """Create a new Nextflow workflow job and enqueue it."""

    payload: Dict[str, Any] = request.model_dump()
    job_id: str = str(uuid.uuid4())

    job_dir: Path = RUN_FILES_PATH / job_id
    nxf_work_dir: Path = job_dir / "_work_dir"
    output_folder: Path = job_dir / "_output"
    job_dir.mkdir(parents=True, exist_ok=True)
    nxf_work_dir.mkdir(parents=True, exist_ok=True)
    output_folder.mkdir(parents=True, exist_ok=True)

    payload["job_dir"] = job_dir
    payload["nxf_work_dir"] = nxf_work_dir
    payload["output_folder"] = output_folder

    job = JobInfo(
        id=job_id,
        status=JobStatus.QUEUED,
        created_at=datetime.now(),
        job_dir=job_dir,
        nxf_work_dir=nxf_work_dir,
        output_folder=output_folder,
        payload=payload,
    )

    with job_store_lock:
        job_store[job_id] = job

    job_queue.put(job_id)

    return job


@app.get(
    "/jobs",
    response_model=List[JobSummary],
)
def list_jobs() -> List[JobSummary]:
    """List all jobs with basic information."""

    with job_store_lock:
        return [
            JobSummary(
                id=j.id,
                status=j.status,
                created_at=j.created_at,
                started_at=j.started_at,
                finished_at=j.finished_at,
            )
            for j in job_store.values()
        ]


@app.get(
    "/jobs/running",
    response_model=List[JobSummary],
)
def list_running_jobs() -> List[JobSummary]:
    """List jobs that are currently running."""

    with job_store_lock:
        return [
            JobSummary(
                id=j.id,
                status=j.status,
                created_at=j.created_at,
                started_at=j.started_at,
                finished_at=j.finished_at,
            )
            for j in job_store.values()
            if j.status == JobStatus.RUNNING
        ]


@app.get(
    "/jobs/completed",
    response_model=List[JobSummary],
)
def list_completed_jobs() -> List[JobSummary]:
    """List jobs that have completed successfully."""

    with job_store_lock:
        return [
            JobSummary(
                id=j.id,
                status=j.status,
                created_at=j.created_at,
                started_at=j.started_at,
                finished_at=j.finished_at,
            )
            for j in job_store.values()
            if j.status == JobStatus.COMPLETED
        ]


@app.get(
    "/jobs/failed",
    response_model=List[JobSummary],
)
def list_failed_jobs() -> List[JobSummary]:
    """List jobs that have failed."""

    with job_store_lock:
        return [
            JobSummary(
                id=j.id,
                status=j.status,
                created_at=j.created_at,
                started_at=j.started_at,
                finished_at=j.finished_at,
            )
            for j in job_store.values()
            if j.status == JobStatus.FAILED
        ]


@app.get(
    "/jobs/{job_id}",
    response_model=JobInfo,
)
def get_job(
    job_id: str,
) -> JobInfo:
    """Retrieve full information for a single job."""

    with job_store_lock:
        job = job_store.get(job_id)

    if job is None:
        raise HTTPException(
            status_code=404,
            detail="Job not found",
        )

    return job


@app.delete(
    "/jobs/{job_id}",
    response_model=JobSummary,
)
def delete_job(
    job_id: str,
) -> JobSummary:
    """Remove a finished job from the in memory store."""

    with job_store_lock:
        job = job_store.get(job_id)
        if job is None:
            raise HTTPException(
                status_code=404,
                detail="Job not found",
            )
        if job.status in {
            JobStatus.QUEUED,
            JobStatus.RUNNING,
        }:
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
        )
        del job_store[job_id]

    return summary


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(
        app=app,
        host="0.0.0.0",
        port=8000,
        reload=False,
    )