#!/usr/bin/env bash
set -e

# ---------- CONFIG ----------
BACKEND_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/backend"
FRONTEND_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/frontend"

BACKEND_HOST="0.0.0.0"
BACKEND_PORT="8000"
# ----------------------------

echo "Starting backend and frontend..."
echo "Backend:  http://localhost:${BACKEND_PORT}"
echo "Frontend: http://localhost:5173"
echo

# Start backend
(
  cd "$BACKEND_DIR"
  echo "[backend] starting..."
  python app.py
) &
BACKEND_PID=$!

# Start frontend
(
  cd "$FRONTEND_DIR"
  echo "[frontend] starting..."
  npm run dev
) &
FRONTEND_PID=$!

# Handle Ctrl+C properly
cleanup() {
  echo
  echo "Stopping backend and frontend..."
  kill "$BACKEND_PID" "$FRONTEND_PID" 2>/dev/null || true
  wait "$BACKEND_PID" "$FRONTEND_PID" 2>/dev/null || true
  echo "Stopped."
}

trap cleanup INT TERM

# Wait for both processes
wait "$BACKEND_PID" "$FRONTEND_PID"
