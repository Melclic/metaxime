import React, { useEffect, useMemo, useState } from "react";
import { useNavigate } from "react-router-dom";

type JobResult = {
  id: string;
  steps: number;
  rp_mean_score: number;
  rp_std_score: number;
};

type JobResultsTableProps = {
  jobId: string;
};

type SortKey = keyof JobResult;
type SortDirection = "asc" | "desc";

export function JobResultsTable({ jobId }: JobResultsTableProps) {
  const [results, setResults] = useState<JobResult[]>([]);
  const [loading, setLoading] = useState<boolean>(true);
  const [error, setError] = useState<string>("");

  const [sortKey, setSortKey] = useState<SortKey>("id");
  const [sortDirection, setSortDirection] = useState<SortDirection>("asc");

  const navigate = useNavigate();

  useEffect(() => {
    if (!jobId) return;

    async function fetchResults() {
      try {
        setLoading(true);
        setError("");

        const response = await fetch(
          `http://localhost:8000/jobs/${encodeURIComponent(jobId)}/results`
        );
        if (!response.ok) {
          throw new Error(`HTTP ${response.status}`);
        }

        const data: JobResult[] = await response.json();
        setResults(data);
      } catch (err: any) {
        console.error(err);
        setError("Could not load job results");
      } finally {
        setLoading(false);
      }
    }

    fetchResults();
  }, [jobId]);

  function handleSort(column: SortKey) {
    if (column === sortKey) {
      setSortDirection(prev => (prev === "asc" ? "desc" : "asc"));
    } else {
      setSortKey(column);
      setSortDirection("asc");
    }
  }

  function handleRowClick(resultId: string) {
    navigate(
      `/jobs/${encodeURIComponent(jobId)}/results/${encodeURIComponent(
        resultId
      )}/pathway`
    );
  }

  const sortedResults = useMemo(() => {
    const copy = [...results];
    copy.sort((a, b) => {
      const aVal = a[sortKey];
      const bVal = b[sortKey];

      if (typeof aVal === "number" && typeof bVal === "number") {
        return sortDirection === "asc" ? aVal - bVal : bVal - aVal;
      }

      const aStr = String(aVal);
      const bStr = String(bVal);
      if (aStr < bStr) return sortDirection === "asc" ? -1 : 1;
      if (aStr > bStr) return sortDirection === "asc" ? 1 : -1;
      return 0;
    });
    return copy;
  }, [results, sortKey, sortDirection]);

  if (!jobId) {
    return <p>Please select a job.</p>;
  }

  if (loading) {
    return <p>Loading results…</p>;
  }

  if (error) {
    return <p style={{ color: "red" }}>{error}</p>;
  }

  if (sortedResults.length === 0) {
    return <p>No results found for this job.</p>;
  }

  return (
    <div style={{ width: "100%", overflowX: "auto" }}>
      <table
        style={{
          borderCollapse: "collapse",
          width: "100%",
          fontSize: "0.9rem",
        }}
      >
        <thead>
          <tr>
            <SortableHeader
              label="Id"
              column="id"
              sortKey={sortKey}
              sortDirection={sortDirection}
              onSort={handleSort}
            />
            <SortableHeader
              label="Steps"
              column="steps"
              sortKey={sortKey}
              sortDirection={sortDirection}
              onSort={handleSort}
            />
            <SortableHeader
              label="RP mean score"
              column="rp_mean_score"
              sortKey={sortKey}
              sortDirection={sortDirection}
              onSort={handleSort}
            />
            <SortableHeader
              label="RP standard score"
              column="rp_std_score"
              sortKey={sortKey}
              sortDirection={sortDirection}
              onSort={handleSort}
            />
          </tr>
        </thead>
        <tbody>
          {sortedResults.map(row => (
            <tr
              key={row.id}
              onClick={() => handleRowClick(row.id)}
              style={rowStyle}
            >
              <td style={cellStyle}>{row.id}</td>
              <td style={cellStyle}>{row.steps}</td>
              <td style={cellStyle}>{formatNumber(row.rp_mean_score)}</td>
              <td style={cellStyle}>{formatNumber(row.rp_std_score)}</td>
            </tr>
          ))}
        </tbody>
      </table>
      <p style={{ fontSize: "0.8rem", opacity: 0.7, marginTop: 8 }}>
        Click a row to view the pathway.
      </p>
    </div>
  );
}

type SortableHeaderProps = {
  label: string;
  column: SortKey;
  sortKey: SortKey;
  sortDirection: SortDirection;
  onSort: (column: SortKey) => void;
};

function SortableHeader({
  label,
  column,
  sortKey,
  sortDirection,
  onSort,
}: SortableHeaderProps) {
  const active = sortKey === column;
  const indicator = active ? (sortDirection === "asc" ? "▲" : "▼") : "";

  return (
    <th
      style={{
        ...cellHeaderStyle,
        cursor: "pointer",
        userSelect: "none",
        whiteSpace: "nowrap",
      }}
      onClick={() => onSort(column)}
    >
      <span style={{ marginRight: 4 }}>{label}</span>
      <span style={{ fontSize: "0.75rem", opacity: active ? 0.8 : 0.3 }}>
        {indicator}
      </span>
    </th>
  );
}

const rowStyle: React.CSSProperties = {
  cursor: "pointer",
};

const cellStyle: React.CSSProperties = {
  border: "1px solid #ddd",
  padding: "6px 8px",
  textAlign: "left",
};

const cellHeaderStyle: React.CSSProperties = {
  ...cellStyle,
  fontWeight: 600,
  backgroundColor: "#f5f5f5",
};

function formatNumber(value: number | null | undefined): string {
  if (value === null || value === undefined) return "";
  return value.toFixed(3);
}