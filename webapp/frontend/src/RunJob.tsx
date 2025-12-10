//import React, { useState, ChangeEvent, DragEvent, FormEvent } from "react";

import React from "react";
import { useState } from "react";
import type { ChangeEvent } from "react";
import type { DragEvent } from "react";
import type { FormEvent } from "react";

import Button from '@mui/material/Button';

type JobResponse = {
  id: string;
  status: string;
  [key: string]: unknown;
};

type UploadModelResponse = {
  model_file: string;
  original_filename: string;
  temporary?: boolean;
};

type UploadRulesResponse = {
  rules_file: string;
  original_filename: string;
  temporary?: boolean;
};

export const RunJob: React.FC = () => {
  const [inchi, setInchi] = useState<string>("");
  const [modelFilePath, setModelFilePath] = useState<string>("");
  const [maxSteps, setMaxSteps] = useState<number | "">(5);
  const [useInchikey2, setUseInchikey2] = useState<boolean>(true);
  const [findAllParentless, setFindAllParentless] = useState<boolean>(true);

  // Advanced options
  const [showAdvanced, setShowAdvanced] = useState<boolean>(false);
  const [rulesFile, setRulesFile] = useState<string>("");
  const [stdMode, setStdMode] = useState<string>("");
  const [topx, setTopx] = useState<number | "">("");
  const [acceptPartialResults, setAcceptPartialResults] =
    useState<boolean>(false);
  const [diameters, setDiameters] = useState<string>("");
  const [ruleType, setRuleType] = useState<string>("");
  const [ramLimit, setRamLimit] = useState<number | "">("");
  const [sourceComp, setSourceComp] = useState<string>("");
  const [targetComp, setTargetComp] = useState<string>("");

  const [selectedFileName, setSelectedFileName] = useState<string>("");
  const [selectedRulesFileName, setSelectedRulesFileName] = useState<string>("");
  const [uploadingModel, setUploadingModel] = useState<boolean>(false);
  const [uploadingRules, setUploadingRules] = useState<boolean>(false);
  const [loading, setLoading] = useState<boolean>(false);
  const [error, setError] = useState<string>("");
  const [jobResponse, setJobResponse] = useState<JobResponse | null>(null);

  // Upload SBML model file
  const uploadModelFile = async (file: File): Promise<void> => {
    setError("");
    setUploadingModel(true);
    try {
      const formData = new FormData();
      formData.append("file", file);

      const res = await fetch("http://localhost:8000/upload_model", {
        method: "POST",
        body: formData,
      });

      if (!res.ok) {
        const text = await res.text();
        throw new Error(`Upload failed HTTP ${res.status}: ${text}`);
      }

      const data = (await res.json()) as UploadModelResponse;
      setSelectedFileName(data.original_filename);
      setModelFilePath(data.model_file);
    } catch (err: any) {
      console.error(err);
      setError(err.message || "Could not upload model file");
    } finally {
      setUploadingModel(false);
    }
  };

  // Upload rules file
  const uploadRulesFile = async (file: File): Promise<void> => {
    setError("");
    setUploadingRules(true);
    try {
      const formData = new FormData();
      formData.append("file", file);

      const res = await fetch("http://localhost:8000/upload_rules", {
        method: "POST",
        body: formData,
      });

      if (!res.ok) {
        const text = await res.text();
        throw new Error(`Rules upload failed HTTP ${res.status}: ${text}`);
      }

      const data = (await res.json()) as UploadRulesResponse;
      setSelectedRulesFileName(data.original_filename);
      // this path will be used as rules_file in the payload
      setRulesFile(data.rules_file);
    } catch (err: any) {
      console.error(err);
      setError(err.message || "Could not upload rules file");
    } finally {
      setUploadingRules(false);
    }
  };

  const handleModelFileInputChange = (
    e: ChangeEvent<HTMLInputElement>,
  ): void => {
    const file = e.target.files?.[0];
    if (!file) return;
    void uploadModelFile(file);
  };

  const handleRulesFileInputChange = (
    e: ChangeEvent<HTMLInputElement>,
  ): void => {
    const file = e.target.files?.[0];
    if (!file) return;
    void uploadRulesFile(file);
  };

  const handleDrop = (e: DragEvent<HTMLDivElement>): void => {
    e.preventDefault();
    e.stopPropagation();

    const file = e.dataTransfer.files?.[0];
    if (!file) return;
    void uploadModelFile(file);
  };

  const handleDragOver = (e: DragEvent<HTMLDivElement>): void => {
    e.preventDefault();
    e.stopPropagation();
  };

  const handleSubmit = async (e: FormEvent): Promise<void> => {
    e.preventDefault();
    setError("");
    setJobResponse(null);

    if (!modelFilePath) {
      setError("Please upload a model file");
      return;
    }
    if (!inchi) {
      setError("Please provide a InChI string.");
      return;
    }

    setLoading(true);

    try {
      const basePayload: Record<string, unknown> = {
        model_file: modelFilePath,
        target_inchi: inchi,
        max_steps: maxSteps === "" ? undefined : maxSteps,
        use_inchikey2: useInchikey2,
        find_all_parentless: findAllParentless,
        // advanced options
        rules_file: rulesFile || undefined,
        std_mode: stdMode || undefined,
        topx: topx === "" ? undefined : topx,
        accept_partial_results: acceptPartialResults || undefined,
        diameters: diameters || undefined,
        rule_type: ruleType || undefined,
        ram_limit: ramLimit === "" ? undefined : ramLimit,
        source_comp: sourceComp || undefined,
        target_comp: targetComp || undefined,
      };

      // Strip undefined / null / empty fields
      const payload = Object.fromEntries(
        Object.entries(basePayload).filter(
          ([, v]) => v !== undefined && v !== null && v !== "",
        ),
      );

      const res = await fetch("http://localhost:8000/jobs", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify(payload),
      });

      if (!res.ok) {
        const text = await res.text();
        throw new Error(
          `Job submit failed HTTP ${res.status}: ${
            text || "unknown error"
          }`,
        );
      }

      const data = (await res.json()) as JobResponse;
      setJobResponse(data);
    } catch (err: any) {
      console.error(err);
      setError(err.message || "Could not submit job");
    } finally {
      setLoading(false);
    }
  };

  return (
    <div
      style={{
        maxWidth: 900,
        margin: "0 auto",
        padding: "16px",
        fontFamily: "Arial, sans-serif",
      }}
    >
      <h1 style={{ marginBottom: 16 }}>Run Job</h1>
        <p>
        The Run Job section allows you to launch a metabolic route prediction workflow.
        Provide the chemical structure of your target molecule and an SBML model
        representing the host organism in which the production pathway should be
        evaluated. Once submitted, the system will compute possible biosynthetic routes
        leading to the desired target within the selected metabolic model.
        </p>
      <form onSubmit={handleSubmit}>
        <div style={{ marginBottom: 12 }}>
          <label
            htmlFor="inchi"
            style={{ display: "block", fontWeight: 600, marginBottom: 4 }}
          >
            Target structure InChI
          </label>
          <textarea
            id="inchi"
            value={inchi}
            onChange={e => setInchi(e.target.value)}
            rows={3}
            style={{
              width: "100%",
              padding: 8,
              fontFamily: "monospace",
              fontSize: 14,
              boxSizing: "border-box",
            }}
            placeholder="Enter a InChI string"
          />
        </div>

        {/* Model upload area */}
        <div style={{ marginBottom: 16 }}>
          <label
            style={{ display: "block", fontWeight: 600, marginBottom: 4 }}
          >
            SBML model file
          </label>

          <div
            onDrop={handleDrop}
            onDragOver={handleDragOver}
            style={{
              border: "2px dashed #ccc",
              borderRadius: 8,
              padding: 16,
              textAlign: "center",
              marginBottom: 8,
              backgroundColor: "#fafafa",
            }}
          >
            <div style={{ marginBottom: 8 }}>
              Drag and drop a model file here
            </div>
            <div style={{ fontSize: 12, opacity: 0.75, marginBottom: 8 }}>
              or click to select a file
            </div>
            <label
              style={{
                display: "inline-block",
                padding: "6px 12px",
                borderRadius: 4,
                backgroundColor: "#1976d2",
                color: "#fff",
                cursor: "pointer",
                fontSize: 13,
              }}
            >
              Choose file
              <input
                type="file"
                accept=".xml,.sbml"
                onChange={handleModelFileInputChange}
                style={{ display: "none" }}
              />
            </label>

            {uploadingModel && (
              <div style={{ marginTop: 8, fontSize: 12 }}>
                Uploading model file…
              </div>
            )}
          </div>

          {selectedFileName && (
            <div style={{ fontSize: 13, marginBottom: 4 }}>
              Uploaded model: <strong>{selectedFileName}</strong>
            </div>
          )}
{/* 
          {modelFilePath && (
            <div style={{ marginTop: 4 }}>
              <label
                htmlFor="model-file-path"
                style={{
                  display: "block",
                  fontWeight: 600,
                  marginBottom: 2,
                  fontSize: 12,
                }}
              >
                Server path used as model_file
              </label>
              <input
                id="model-file-path"
                type="text"
                value={modelFilePath}
                onChange={e => setModelFilePath(e.target.value)}
                style={{
                  width: "100%",
                  padding: 6,
                  boxSizing: "border-box",
                  fontFamily: "monospace",
                  fontSize: 12,
                }}
              />
            </div>
          )}
           */}
        </div>

        {/* Max steps */}
        <div style={{ marginBottom: 12 }}>
          <label
            htmlFor="max-steps"
            style={{ display: "block", fontWeight: 600, marginBottom: 4 }}
          >
            Maximum steps
          </label>
          <input
            id="max-steps"
            type="number"
            min={1}
            value={maxSteps}
            onChange={e =>
              setMaxSteps(e.target.value === "" ? "" : Number(e.target.value))
            }
            style={{
              width: 150,
              padding: 6,
              boxSizing: "border-box",
            }}
            placeholder="e.g. 5"
          />
        </div>

        {/* Basic flags */}
        <div style={{ marginBottom: 12 }}>
          <label style={{ display: "block", fontWeight: 600, marginBottom: 4 }}>
            Basic options
          </label>
          <label
            style={{
              display: "inline-flex",
              alignItems: "center",
              marginRight: 16,
            }}
          >
            <input
              type="checkbox"
              checked={useInchikey2}
              onChange={e => setUseInchikey2(e.target.checked)}
              style={{ marginRight: 6 }}
            />
            Use InChIKey2 fallback matching
          </label>

          <label
            style={{
              display: "inline-flex",
              alignItems: "center",
            }}
          >
            <input
              type="checkbox"
              checked={findAllParentless}
              onChange={e => setFindAllParentless(e.target.checked)}
              style={{ marginRight: 6 }}
            />
            Find all parentless metabolites
          </label>
        </div>

        {/* Advanced options (collapsible) */}
        <div
          style={{
            border: "1px solid #ddd",
            borderRadius: 8,
            padding: 10,
            marginBottom: 16,
            backgroundColor: "#fafafa",
          }}
        >
          <Button
            variant="text"
            color="primary"
            size="small"
            
            onClick={() => setShowAdvanced(v => !v)}

            // type="button"
            // style={{
            //   border: "none",
            //   background: "none",
            //   padding: 0,
            //   margin: 0,
            //   fontSize: 14,
            //   fontWeight: 600,
            //   cursor: "pointer",
            //   display: "flex",
            //   alignItems: "center",
            //   marginBottom: showAdvanced ? 8 : 0,
            // }}
          >
            <span style={{ marginRight: 6 }}>
              {showAdvanced ? "▼" : "▶"}
            </span>
            Advanced options
          </Button>

          {showAdvanced && (
            <div style={{ marginTop: 4 }}>
              {/* Rules file with simple upload button */}
              <div style={{ marginBottom: 8 }}>
                <label
                  style={{
                    display: "block",
                    fontWeight: 600,
                    marginBottom: 2,
                  }}
                >
                  Rules file
                </label>
                {/* <input
                  type="text"
                  value={rulesFile}
                  onChange={e => setRulesFile(e.target.value)}
                  style={{
                    width: "100%",
                    padding: 6,
                    boxSizing: "border-box",
                    marginBottom: 4,
                  }}
                  placeholder="e.g. /path/to/rules.tsv"
                /> */}
                <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
                  <label
                    style={{
                      display: "inline-block",
                      padding: "4px 10px",
                      borderRadius: 4,
                      backgroundColor: "#1976d2",
                      color: "#fff",
                      cursor: "pointer",
                      fontSize: 12,
                    }}
                  >
                    Choose file
                    <input
                      type="file"
                      onChange={handleRulesFileInputChange}
                      style={{ display: "none" }}
                    />
                  </label>
                  {uploadingRules && (
                    <span style={{ fontSize: 12 }}>Uploading rules…</span>
                  )}
                </div>
                {selectedRulesFileName && (
                  <div style={{ fontSize: 12, marginTop: 4 }}>
                    Uploaded rules:{" "}
                    <strong>{selectedRulesFileName}</strong>
                  </div>
                )}
              </div>

              <div style={{ marginBottom: 8 }}>
                <label
                  style={{
                    display: "block",
                    fontWeight: 600,
                    marginBottom: 2,
                  }}
                >
                  Standardisation mode
                </label>
                <input
                  type="text"
                  value={stdMode}
                  onChange={e => setStdMode(e.target.value)}
                  style={{
                    width: "100%",
                    padding: 6,
                    boxSizing: "border-box",
                  }}
                  placeholder="e.g. H added + Aromatized"
                />
              </div>

              <div style={{ marginBottom: 8 }}>
                <label
                  style={{
                    display: "block",
                    fontWeight: 600,
                    marginBottom: 2,
                  }}
                >
                  Top X pathways
                </label>
                <input
                  type="number"
                  min={1}
                  value={topx}
                  onChange={e =>
                    setTopx(
                      e.target.value === "" ? "" : Number(e.target.value),
                    )
                  }
                  style={{
                    width: 150,
                    padding: 6,
                    boxSizing: "border-box",
                  }}
                  placeholder="e.g. 10"
                />
              </div>

              <div style={{ marginBottom: 8 }}>
                <label
                  style={{
                    display: "block",
                    fontWeight: 600,
                    marginBottom: 2,
                  }}
                >
                  Diameters
                </label>
                <input
                  type="text"
                  value={diameters}
                  onChange={e => setDiameters(e.target.value)}
                  style={{
                    width: "100%",
                    padding: 6,
                    boxSizing: "border-box",
                  }}
                  placeholder="e.g. 2,4,6"
                />
              </div>

              <div style={{ marginBottom: 8 }}>
                <label
                  style={{
                    display: "block",
                    fontWeight: 600,
                    marginBottom: 2,
                  }}
                >
                  Rule type
                </label>
                <input
                  type="text"
                  value={ruleType}
                  onChange={e => setRuleType(e.target.value)}
                  style={{
                    width: "100%",
                    padding: 6,
                    boxSizing: "border-box",
                  }}
                  placeholder="e.g. all, forward, reverse"
                />
              </div>

              <div style={{ marginBottom: 8 }}>
                <label
                  style={{
                    display: "block",
                    fontWeight: 600,
                    marginBottom: 2,
                  }}
                >
                  RAM limit in GB
                </label>
                <input
                  type="number"
                  min={1}
                  value={ramLimit}
                  onChange={e =>
                    setRamLimit(
                      e.target.value === "" ? "" : Number(e.target.value),
                    )
                  }
                  style={{
                    width: 150,
                    padding: 6,
                    boxSizing: "border-box",
                  }}
                  placeholder="e.g. 32"
                />
              </div>

              <div style={{ marginBottom: 8 }}>
                <label
                  style={{
                    display: "block",
                    fontWeight: 600,
                    marginBottom: 2,
                  }}
                >
                  Source compartment
                </label>
                <input
                  type="text"
                  value={sourceComp}
                  onChange={e => setSourceComp(e.target.value)}
                  style={{
                    width: "100%",
                    padding: 6,
                    boxSizing: "border-box",
                  }}
                  placeholder="e.g. c"
                />
              </div>

              <div style={{ marginBottom: 8 }}>
                <label
                  style={{
                    display: "block",
                    fontWeight: 600,
                    marginBottom: 2,
                  }}
                >
                  Target compartment
                </label>
                <input
                  type="text"
                  value={targetComp}
                  onChange={e => setTargetComp(e.target.value)}
                  style={{
                    width: "100%",
                    padding: 6,
                    boxSizing: "border-box",
                  }}
                  placeholder="e.g. c"
                />
              </div>

              <div style={{ marginBottom: 8 }}>
                <label
                  style={{
                    display: "block",
                    fontWeight: 600,
                    marginBottom: 2,
                  }}
                >
                  Advanced flags
                </label>
                <label
                  style={{
                    display: "flex",
                    alignItems: "center",
                    marginBottom: 4,
                  }}
                >
                  <input
                    type="checkbox"
                    checked={acceptPartialResults}
                    onChange={e =>
                      setAcceptPartialResults(e.target.checked)
                    }
                    style={{ marginRight: 6 }}
                  />
                  Accept partial results
                </label>
              </div>
            </div>
          )}
        </div>

        {/* Error + response */}
        {error && (
          <div
            style={{
              marginBottom: 12,
              color: "red",
              fontSize: 13,
              whiteSpace: "pre-wrap",
            }}
          >
            {error}
          </div>
        )}

        {jobResponse && (
          <div
            style={{
              marginBottom: 12,
              fontSize: 13,
              padding: 8,
              borderRadius: 6,
              border: "1px solid #d0e0d0",
              backgroundColor: "#f3fff3",
            }}
          >
            <div>
              <strong>Job created</strong>
            </div>
            <div>Id: {jobResponse.id}</div>
            <div>Status: {String(jobResponse.status)}</div>
          </div>
        )}

        <button
          disabled={loading || uploadingModel || uploadingRules}
          //variant="contained"
          type="submit"
          style={{
            padding: "8px 16px",
            fontSize: 14,
            borderRadius: 6,
            border: "1px solid #1976d2",
            backgroundColor:
              loading || uploadingModel || uploadingRules
                ? "#90caf9"
                : "#1976d2",
            color: "#fff",
            cursor:
              loading || uploadingModel || uploadingRules
                ? "default"
                : "pointer",
          }}
        >
          {loading ? "Submitting…" : "Run job"}
        </button>
      </form>
    </div>
  );
};