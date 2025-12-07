import { useEffect, useState } from "react";

type JobSummary = {
  id: string;
  status: string;
  created_at: string;
  started_at: string | null;
  finished_at: string | null;
  target_inchi: string | null;
};

export function JobsTable() {
  const [jobs, setJobs] = useState<JobSummary[]>([]);
  const [loading, setLoading] = useState<boolean>(true);
  const [error, setError] = useState<string>("");

  useEffect(() => {
    async function fetchJobs() {
      try {
        setLoading(true);
        setError("");
        const response = await fetch("http://localhost:8000/jobs");
        if (!response.ok) {
          throw new Error(`HTTP ${response.status}`);
        }
        const data: JobSummary[] = await response.json();
        setJobs(data);
      } catch (err: any) {
        console.error(err);
        setError("Could not load jobs");
      } finally {
        setLoading(false);
      }
    }

    fetchJobs();
  }, []);

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
    <table style={{ borderCollapse: "collapse", width: "100%" }}>
      <thead>
        <tr>
          <th style={cellStyle}>Id</th>
          <th style={cellStyle}>Status</th>
          <th style={cellStyle}>Created</th>
          <th style={cellStyle}>Started</th>
          <th style={cellStyle}>Finished</th>
          <th style={cellStyle}>InChI</th>
        </tr>
      </thead>
      <tbody>
        {jobs.map(job => (
          <tr key={job.id}>
            <td style={cellStyle}>{job.id}</td>
            <td style={cellStyle}>{job.status}</td>
            <td style={cellStyle}>{job.created_at}</td>
            <td style={cellStyle}>{job.started_at}</td>
            <td style={cellStyle}>{job.finished_at}</td>
            <td style={cellStyle}>{job.target_inchi}</td>
          </tr>
        ))}
      </tbody>
    </table>
  );
}

const cellStyle: React.CSSProperties = {
  border: "1px solid #ccc",
  padding: "4px 8px",
  textAlign: "left"
};