// App.tsx

import { BrowserRouter, Routes, Route, useParams } from "react-router-dom";



import { JobsTiles } from "./JobsTiles";
import { JobResultsTable } from "./JobResultsTable";
import { PathwayViewer } from "./PathwayViewer";
import { RunJob } from "./RunJob";

import CssBaseline from '@mui/material/CssBaseline';
//import Divider from '@mui/material/Divider';

import AppTheme from './shared-theme/AppTheme';
import AppAppBar from './components/AppAppBar';
import { Toolbar } from "@mui/material";
import Box from '@mui/material/Box';



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
    <AppTheme>
      <CssBaseline enableColorScheme />

      <Box sx={{ pt: 3 }}>
        <Toolbar /> {/* pushes content below AppBar */}

        <Box
          sx={{
            px: 4,            // left + right padding
            width: "100%",
            boxSizing: "border-box",
          }}
        >

        <BrowserRouter>
          <AppAppBar />
          <Routes>
            <Route path="/" element={<JobsTiles />} />
            <Route path="/jobs/:jobId/results" element={<JobResultsRouteWrapper />} />
            <Route path="/jobs/:jobId/results/:resultId" element={<PathwayViewerRouteWrapper />}
            />
            <Route path="/run" element={<RunJob />} />
          </Routes>
        </BrowserRouter>
      
      </Box>
      </Box>

      </AppTheme>
  );
}

export default App;