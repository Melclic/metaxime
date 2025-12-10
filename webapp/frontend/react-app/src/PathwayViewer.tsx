// PathwayViewer.tsx
import React, { useEffect, useState, useMemo } from "react";
import { PathwayNetwork } from "./PathwayNetwork";
import type { PathwayGraph } from "./PathwayNetwork";

type PathwayViewerProps = {
  jobId: string;
  resultId: string;
};

type ReactionNode = {
  id: string;
  type: "reaction";
  annotation?: {
    [key: string]: unknown;
  };
};

export const PathwayViewer: React.FC<PathwayViewerProps> = ({
  jobId,
  resultId,
}) => {
  const [graph, setGraph] = useState<PathwayGraph | null>(null);
  const [loading, setLoading] = useState<boolean>(true);
  const [error, setError] = useState<string>("");

  useEffect(() => {
    async function fetchGraph() {
      try {
        setLoading(true);
        setError("");

        const res = await fetch(
          `http://localhost:8000/jobs/${encodeURIComponent(
            jobId,
          )}/results/${encodeURIComponent(resultId)}`,
        );
        if (!res.ok) {
          throw new Error(`HTTP ${res.status}`);
        }

        const data = (await res.json()) as PathwayGraph;
        setGraph(data);
      } catch (err) {
        console.error(err);
        setError("Could not load pathway");
      } finally {
        setLoading(false);
      }
    }

    fetchGraph();
  }, [jobId, resultId]);

  const reactions: ReactionNode[] = useMemo(() => {
    if (!graph) return [];
    return graph.nodes.filter(
      (n: any) => n.type === "reaction",
    ) as ReactionNode[];
  }, [graph]);

  if (loading) {
    return <p>Loading pathway…</p>;
  }

  if (error) {
    return <p style={{ color: "red" }}>{error}</p>;
  }

  if (!graph) {
    return <p>No pathway data found.</p>;
  }

  return (
    <div
        style={{
        width: "100%",
        marginLeft: "0",
        boxSizing: "border-box",
        padding: 16,
        }}
    >   

    <h2>
    Pathway for job {jobId}, result {resultId}
    </h2>

      <h3 style={{ marginTop: 24, marginBottom: 8 }}>Reactions</h3>
      {reactions.length === 0 ? (
        <p>No reactions found.</p>
      ) : (
        <div style={{ width: "100%" }}>
          <table
            style={{
              borderCollapse: "collapse",
              width: "100%",
              fontSize: "0.9rem",
              fontFamily: "Arial, sans-serif",
            }}
          >
            <thead>
              <tr>
                <th style={thStyle}>ID</th>
                <th style={thStyle}>EC code</th>
                <th style={thStyle}>RP ID</th>
                <th style={thStyle}>RP step</th>
                <th style={thStyle}>RP score</th>
              </tr>
            </thead>
            <tbody>
              {reactions.map(r => {
                const ann = r.annotation ?? {};
                const ecCodes = (ann["ec-code"] as string[]) ?? [];
                const rpId = ann["rp_id"] as string | undefined;
                const rpStep = ann["rp_step"] as number | undefined;
                const rpScore = ann["rp_score"] as number | undefined;

                return (
                  <tr key={r.id}>
                    <td style={tdStyle}>{r.id}</td>
                    <td style={tdStyle}>{ecCodes.join(", ")}</td>
                    <td style={tdStyle}>{rpId ?? "–"}</td>
                    <td style={tdStyle}>
                      {rpStep !== undefined ? rpStep : "–"}
                    </td>
                    <td style={tdStyle}>
                      {rpScore !== undefined
                        ? rpScore.toFixed(3)
                        : "–"}
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
      )}

    <div style={{ width: "100%", marginTop: 16, marginBottom: 16 }}>
        <PathwayNetwork graph={graph} />
    </div>

    </div>
  );
};

const thStyle: React.CSSProperties = {
  border: "1px solid #ddd",
  padding: "6px 8px",
  textAlign: "left",
  backgroundColor: "#f5f5f5",
  fontWeight: 600,
};

const tdStyle: React.CSSProperties = {
  border: "1px solid #ddd",
  padding: "6px 8px",
  textAlign: "left",
};