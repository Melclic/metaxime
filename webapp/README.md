# Metaxime App

**Metaxime** is a full‑stack application for the prediction of metabolic routes enabling the production of a **molecule of interest** within a **strain (metabolic model) of interest**.

The application provides:
- A **FastAPI backend** to run metabolic pathway prediction jobs
- A **React frontend** to submit jobs, monitor status, and visualize pathways
- A **Nextflow** workflow executed behind the scenes to perform the metabolic route prediction

Docker **must be installed** on your system.  
While Nextflow executes locally, several dependencies require Docker to exist on the system even if not actively used to containerize the backend.

---

## Features

- Submit new computational pathway prediction jobs  
- Upload SBML metabolic models  
- Define target molecules with SMILES or InChI  
- Configure optional advanced parameters  
- Monitor job execution  
- View metabolic pathways interactively  
- Export pathway diagrams to SVG  

---

## Requirements

- **Docker** (must be installed)
- **Conda** (Anaconda or Miniconda recommended)
- **Nextflow**
- macOS or Linux recommended

---

## Installation

### 1. Create and activate the Conda environment

```bash
conda env create -f environment.yml
conda activate metaxime_app
```

### 2. Install frontend dependencies

```bash
cd frontend
npm install
cd ..
```

---

## Running the Application

Use the provided script to launch both backend and frontend:

```bash
./run_app.sh
```

This starts:
- **FastAPI backend** at `http://localhost:8000`
- **Vite/React frontend** at `http://localhost:5173`

Open the UI:

```
http://localhost:5173
```

---

## Project Structure

```
metaxime-app/
│
├── backend/
│   ├── app.py               # FastAPI backend
│   ├── jobs/                # Job execution & Nextflow integration
│   ├── workflows/           # Nextflow workflow handling
│   └── utils/
│
├── frontend/
│   ├── src/
│   │   ├── pages/           # React UI pages
│   │   ├── components/      # Shared components
│   │   └── App.tsx          # Frontend routing
│   └── package.json
│
├── nextflow/
│   ├── main.nf              # Nextflow pipeline
│   └── nextflow.config
│
├── environment.yml          # Conda env definition
├── run_app.sh               # Launch script
└── README.md
```

---

## How Metaxime Works

1. The user submits a job via the UI, providing:
   - A target molecule
   - An SBML model
   - Optional parameters such as step count, etc...

2. The **FastAPI backend**:
   - Creates a temporary job directory
   - Constructs a valid Nextflow command
   - Executes the workflow locally

3. **Nextflow** performs:
   - Sink generation
   - RetroPath2
   - Monocomponent reactions to full reactions

4. The backend exposes results as JSON.

5. The frontend presents:
   - Job list and status
   - Pathway result tables
   - Interactive metabolic route visualizations

---

## Notes

- Docker must be available on the system.
- The backend itself is not containerized; Nextflow runs natively.
- Job results and logs persist until removed manually.

---

## License