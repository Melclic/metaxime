// App.tsx
import React from "react";
import { BrowserRouter, Routes, Route, useParams } from "react-router-dom";
import { JobsTiles } from "./JobsTiles";
import { JobResultsTable } from "./JobResultsTable";
import { PathwayViewer } from "./PathwayViewer";


function JobResultsRouteWrapper() {
  const { jobId } = useParams<{ jobId: string }>();
  if (!jobId) {
    return <p>Invalid job id</p>;
  }
  return <JobResultsTable jobId={jobId} />;
}

function PathwayViewerRouteWrapper() {
  const { jobId, resultId } = useParams<{ jobId: string; resultId: string }>();

  if (!jobId || !resultId) {
    return <p>Invalid job or result id</p>;
  }

  return <PathwayViewer jobId={jobId} resultId={resultId} />;
}

export function App() {
  return (
    <BrowserRouter>
      <Routes>
        <Route path="/" element={<JobsTiles />} />
        <Route path="/jobs/:jobId/results" element={<JobResultsRouteWrapper />} />
        <Route
          path="/jobs/:jobId/results/:resultId/pathway"
          element={<PathwayViewerRouteWrapper />}
        />
      </Routes>
    </BrowserRouter>
  );
}

export default App;