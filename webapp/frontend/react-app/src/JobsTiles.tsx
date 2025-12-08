// JobsTiles.tsx
import React, { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";
// @ts-ignore
import { SmilesSVG } from "./SmilesSVG";

type JobStatusResponse = {
  job_id: string;
  status: string;
  rp2_status_code: number | null;
};

const RP2_STATUS_MESSAGES: Record<number, string> = {
  0: "OK",
  10: "Timeout error",
  20: "Memory error",
  30: "Source in sink",
  31: "Source in sink not found",
  40: "No result produced",
  50: "Operating system error",
  60: "Out of RAM"
};

type JobSummary = {
  id: string;
  status: string;
  created_at: string;
  started_at: string | null;
  finished_at: string | null;
  // optional, in case backend does not always provide it
  target_smiles?: string | null;
};

type FieldProps = {
  label: string;
  value: string;
  multiline?: boolean;
};

function useResponsiveColumns(): number {
  const [columns, setColumns] = useState<number>(1);

  useEffect(() => {
    function updateColumns() {
      const width = window.innerWidth;

      if (width < 640) {
        setColumns(1);
      } else if (width < 1024) {
        setColumns(2);
      } else if (width < 1440) {
        setColumns(3);
      } else {
        setColumns(4);
      }
    }

    updateColumns();
    window.addEventListener("resize", updateColumns);
    return () => window.removeEventListener("resize", updateColumns);
  }, []);

  return columns;
}

export function JobsTiles() {
  const [jobs, setJobs] = useState<JobSummary[]>([]);
  const [jobErrors, setJobErrors] = useState<Record<string, string>>({});
  const [loading, setLoading] = useState<boolean>(true);

  const [error, setError] = useState<string>("");

  const columns = useResponsiveColumns();
  const navigate = useNavigate();

  useEffect(() => {
    async function fetchJobs() {
      try {
        //setLoading(true);
        //setError("");
        const response = await fetch("http://localhost:8000/jobs");
        if (!response.ok) {
          throw new Error(`HTTP ${response.status}`);
        }
        const data: JobSummary[] = await response.json();
        setJobs(data);
      } catch (err: any) {
        console.error(err);
        setError("Could not load jobs");
      } //finally {
      //   setLoading(false);
      // }
    }

    fetchJobs();
  }, []);

  useEffect(() => {
    async function fetchFailedJobStatuses() {
      const failedJobs = jobs.filter(j => j.status === "failed");

      if (failedJobs.length === 0) {
        setLoading(false);
        return;
      }

      const entries = await Promise.all(
        failedJobs.map(async job => {
          try {
            const res = await fetch(
              `http://localhost:8000/jobs/${encodeURIComponent(
                job.id,
              )}/status`,
            );
            if (!res.ok) {
              throw new Error(`HTTP ${res.status}`);
            }
            const data: JobStatusResponse = await res.json();

            if (data.rp2_status_code !== null) {
              const msg =
                RP2_STATUS_MESSAGES[data.rp2_status_code] ??
                `Unknown error (${data.rp2_status_code})`;
              return [job.id, msg] as const;
            }
            return [job.id, "Unknown failure"] as const;
          } catch (err) {
            console.error(err);
            return [job.id, "Unknown failure"] as const;
          }
        }),
      );

      setJobErrors(Object.fromEntries(entries));
      setLoading(false);
    }

    if (jobs.length > 0) {
      fetchFailedJobStatuses();
    } else {
      setLoading(false);
    }
  }, [jobs]);

  if (loading) {
    return <p>Loading jobs…</p>;
  }

  if (error) {
    return <p style={{ color: "red" }}>{error}</p>;
  }

  if (jobs.length === 0) {
    return <p>No jobs found</p>;
  }

  return (
    <div
      style={{
        width: "100%",
        boxSizing: "border-box",
        padding: "16px 0",
      }}
    >
      <div
        style={{
          display: "grid",
          gridTemplateColumns: `repeat(${columns}, minmax(0, 1fr))`,
          gap: "16px",
          alignItems: "stretch",
        }}
      >
        {jobs.map(job => {
          const isFailed = job.status === "failed";
          const isCompleted = job.status === "completed";
          const errorMessage = jobErrors[job.id];

          return (
            <div
              key={job.id}
              style={{
                ...cardStyle,
                cursor: isCompleted ? "pointer" : "default",
                opacity: isFailed ? 0.7 : 1,
              }}
              onClick={() => {
                if (isCompleted) {
                  navigate(
                    `/jobs/${encodeURIComponent(job.id)}/results`
                  );
                }
              }}
              role={isCompleted ? "button" : undefined}
            >
              <div style={cardHeaderStyle}>
                <span style={cardTitleStyle}>Job</span>
                <span style={statusBadgeStyle(job.status)}>
                  {job.status}
                </span>
              </div>

              <div style={cardBodyStyle}>
                <Field label="Id" value={job.id} />
                <Field label="Created" value={formatDate(job.created_at)} />
                <Field label="Started" value={formatDate(job.started_at)} />
                <Field label="Finished" value={formatDate(job.finished_at)} />

                {isFailed && errorMessage && (
                  <Field
                    label="Error"
                    value={errorMessage}
                    multiline
                  />
                )}
              </div>

              <div style={structureContainerStyle}>
                <div
                  style={{
                    fontSize: "0.8rem",
                    opacity: 0.7,
                    marginBottom: 4,
                  }}
                >
                  Structure
                </div>
                {job.target_smiles ? (
                  <SmilesSVG
                    smiles={job.target_smiles}
                    width={260}
                    height={160}
                  />
                ) : (
                  <div style={{ fontSize: "0.9rem", opacity: 0.7 }}>
                    No structure available
                  </div>
                )}
              </div>
            </div>
          );
        })}
      </div>
    </div>
  );

}

function Field({ label, value, multiline }: FieldProps) {
  return (
    <div style={{ marginBottom: "6px" }}>
      <div style={{ fontSize: "0.8rem", opacity: 0.7 }}>{label}</div>
      <div
        style={{
          fontSize: "0.9rem",
          whiteSpace: multiline ? "normal" : "nowrap",
          overflow: "hidden",
          textOverflow: "ellipsis",
        }}
        title={value}
      >
        {value}
      </div>
    </div>
  );
}

function formatDate(value: string | null): string {
  if (!value) {
    return "Not available";
  }
  const date = new Date(value);
  if (Number.isNaN(date.getTime())) {
    return value;
  }
  return date.toLocaleString();
}

const cardStyle: React.CSSProperties = {
  borderRadius: 12,
  border: "1px solid #e0e0e0",
  padding: "12px 14px",
  boxShadow: "0 2px 4px rgba(0,0,0,0.04)",
  backgroundColor: "#ffffff",
  display: "flex",
  flexDirection: "column",
  justifyContent: "space-between",
  cursor: "pointer",
  transition: "box-shadow 0.15s ease, transform 0.15s ease",
};

const cardHeaderStyle: React.CSSProperties = {
  display: "flex",
  alignItems: "center",
  justifyContent: "space-between",
  marginBottom: 8,
};

const cardTitleStyle: React.CSSProperties = {
  fontWeight: 600,
  fontSize: "0.95rem",
};

const cardBodyStyle: React.CSSProperties = {
  display: "flex",
  flexDirection: "column",
};

const structureContainerStyle: React.CSSProperties = {
  marginTop: 8,
};

const statusBadgeStyle = (status: string): React.CSSProperties => {
  let background = "#e0e0e0";
  let color = "#333";
  const s = status.toLowerCase();

  if (s === "completed") {
    background = "#e0f7e9";
    color = "#137333";
  } else if (s === "running") {
    background = "#e3f2fd";
    color = "#0b5394";
  } else if (s === "failed") {
    background = "#ffebee";
    color = "#b71c1c";
  } else if (s === "queued") {
    background = "#fff8e1";
    color = "#8d6e00";
  }

  return {
    borderRadius: 999,
    padding: "2px 10px",
    fontSize: "0.75rem",
    fontWeight: 600,
    backgroundColor: background,
    color,
    textTransform: "uppercase",
  };
};