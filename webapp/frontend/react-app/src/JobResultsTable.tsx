// JobResultsTable.tsx
import React, { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";
import {
  Box,
  CircularProgress,
  Paper,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  TableSortLabel,
  Typography,
} from "@mui/material";

type JobResultsTableProps = {
  jobId: string;
};

type JobResult = {
  id: string;
  steps: number;
  rp_mean_score?: number;
  rp_std_score?: number;
};

type Order = "asc" | "desc";

// Comparators specialised for JobResult so TS is happy with optional fields
function descendingComparator(
  a: JobResult,
  b: JobResult,
  orderBy: keyof JobResult,
): number {
  const av = a[orderBy];
  const bv = b[orderBy];

  if (av == null && bv == null) return 0;
  if (av == null) return 1;
  if (bv == null) return -1;

  if (bv < av) return -1;
  if (bv > av) return 1;
  return 0;
}

function getComparator(
  order: Order,
  orderBy: keyof JobResult,
): (a: JobResult, b: JobResult) => number {
  return order === "desc"
    ? (a, b) => descendingComparator(a, b, orderBy)
    : (a, b) => -descendingComparator(a, b, orderBy);
}

function stableSort(
  array: JobResult[],
  comparator: (a: JobResult, b: JobResult) => number,
): JobResult[] {
  const stabilized: Array<[JobResult, number]> = array.map(
    (el, index) => [el, index],
  );
  stabilized.sort((a, b) => {
    const order = comparator(a[0], b[0]);
    if (order !== 0) return order;
    return a[1] - b[1];
  });
  return stabilized.map(el => el[0]);
}

export const JobResultsTable: React.FC<JobResultsTableProps> = ({
  jobId,
}) => {
  const navigate = useNavigate();

  const [results, setResults] = useState<JobResult[]>([]);
  const [loading, setLoading] = useState<boolean>(true);
  const [error, setError] = useState<string>("");

  const [order, setOrder] = useState<Order>("asc");
  const [orderBy, setOrderBy] = useState<keyof JobResult>("id");

  useEffect(() => {
    async function fetchResults() {
      try {
        setLoading(true);
        setError("");

        const res = await fetch(
          `http://localhost:8000/jobs/${encodeURIComponent(jobId)}/results`,
        );
        if (!res.ok) {
          throw new Error(`HTTP ${res.status}`);
        }
        const data: JobResult[] = await res.json();
        setResults(data);
      } catch (err: any) {
        console.error(err);
        setError("Could not load results");
      } finally {
        setLoading(false);
      }
    }

    void fetchResults();
  }, [jobId]);

  const handleRequestSort = (property: keyof JobResult) => {
    const isAsc = orderBy === property && order === "asc";
    setOrder(isAsc ? "desc" : "asc");
    setOrderBy(property);
  };

  const handleRowClick = (resultId: string) => {
    navigate(
      `/jobs/${encodeURIComponent(
        jobId,
      )}/results/${encodeURIComponent(resultId)}`,
    );
  };

  if (loading) {
    return (
      <Box
        sx={{
          mt: 4,
          display: "flex",
          justifyContent: "center",
        }}
      >
        <CircularProgress />
      </Box>
    );
  }

  if (error) {
    return (
      <Box sx={{ mt: 4 }}>
        <Typography color="error">{error}</Typography>
      </Box>
    );
  }

  if (results.length === 0) {
    return (
      <Box sx={{ mt: 4 }}>
        <Typography>No results found for this job.</Typography>
      </Box>
    );
  }

  return (
    <Box sx={{ mt: 2 }}>
      <Typography variant="h5" gutterBottom>
        Pathway results
      </Typography>
      <Typography variant="body2" sx={{ mb: 2 }}>
        Click on a row to open the pathway network for that result.
      </Typography>

      <TableContainer
        component={Paper}
        sx={{
          maxHeight: "70vh",
        }}
      >
        <Table stickyHeader size="small" aria-label="job results table">
          <TableHead>
            <TableRow>
              <TableCell sortDirection={orderBy === "id" ? order : false}>
                <TableSortLabel
                  active={orderBy === "id"}
                  direction={orderBy === "id" ? order : "asc"}
                  onClick={() => handleRequestSort("id")}
                >
                  Id
                </TableSortLabel>
              </TableCell>

              <TableCell
                sortDirection={orderBy === "steps" ? order : false}
                align="right"
              >
                <TableSortLabel
                  active={orderBy === "steps"}
                  direction={orderBy === "steps" ? order : "asc"}
                  onClick={() => handleRequestSort("steps")}
                >
                  Steps
                </TableSortLabel>
              </TableCell>

              <TableCell
                sortDirection={orderBy === "rp_mean_score" ? order : false}
                align="right"
              >
                <TableSortLabel
                  active={orderBy === "rp_mean_score"}
                  direction={orderBy === "rp_mean_score" ? order : "asc"}
                  onClick={() => handleRequestSort("rp_mean_score")}
                >
                  Mean score
                </TableSortLabel>
              </TableCell>

              <TableCell
                sortDirection={orderBy === "rp_std_score" ? order : false}
                align="right"
              >
                <TableSortLabel
                  active={orderBy === "rp_std_score"}
                  direction={orderBy === "rp_std_score" ? order : "asc"}
                  onClick={() => handleRequestSort("rp_std_score")}
                >
                  Std score
                </TableSortLabel>
              </TableCell>
            </TableRow>
          </TableHead>

          <TableBody>
            {stableSort(results, getComparator(order, orderBy)).map(
              row => (
                <TableRow
                  hover
                  key={row.id}
                  onClick={() => handleRowClick(row.id)}
                  sx={{
                    cursor: "pointer",
                  }}
                >
                  <TableCell component="th" scope="row">
                    {row.id}
                  </TableCell>
                  <TableCell align="right">{row.steps}</TableCell>
                  <TableCell align="right">
                    {row.rp_mean_score != null
                      ? row.rp_mean_score.toFixed(3)
                      : "–"}
                  </TableCell>
                  <TableCell align="right">
                    {row.rp_std_score != null
                      ? row.rp_std_score.toFixed(3)
                      : "–"}
                  </TableCell>
                </TableRow>
              ),
            )}
          </TableBody>
        </Table>
      </TableContainer>
    </Box>
  );
};