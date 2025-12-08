import { BrowserRouter, Routes, Route, useParams } from "react-router-dom";
import { JobsTiles } from "./JobsTiles";
import { JobResultsTable } from "./JobResultsTable";


function JobResultsTableRouteWrapper() {
  const { jobId } = useParams<{ jobId: string }>();

  if (!jobId) {
    return <p>Invalid job id</p>;
  }

  return <JobResultsTable jobId={jobId} />;
}

function App() {
  return (
    <div style={{ padding: 24 }}>
      <h1>Jobs</h1>
        <BrowserRouter>
          <Routes>
            <Route path="/" element={<JobsTiles />} />
            <Route path="/jobs/:jobId/results" element={<JobResultsTableRouteWrapper />} />
          </Routes>
        </BrowserRouter>
    </div>
  );
}

export default App;