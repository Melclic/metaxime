// PathwayNetwork.tsx
import React, { useEffect, useRef, useState, useMemo } from "react";
import { createRoot } from "react-dom/client";
import * as d3 from "d3";
import * as dagreD3 from "dagre-d3-es";
import { SmilesSVG } from "./SmilesSVG";

const chemWidth = 140;
const chemHeight = 100;

type MetaboliteNode = {
  id: string;
  type: "metabolite";
  name?: string;
  topology: string;
  is_cofactor?: boolean;
  annotation?: {
    smiles?: string;
    [key: string]: unknown;
  };
  [key: string]: unknown;
};

type ReactionNode = {
  id: string;
  type: "reaction";
  name?: string;
  topology: string;
  annotation?: {
    [key: string]: unknown;
  };
  [key: string]: unknown;
};

type PathwayNode = MetaboliteNode | ReactionNode;

type PathwayLink = {
  source: string;
  target: string;
  role: "reactant" | "product";
  stoichiometry: number;
  [key: string]: unknown;
};

export type PathwayGraph = {
  directed: boolean;
  multigraph: boolean;
  nodes: PathwayNode[];
  links: PathwayLink[];
  id?: string;
  steps?: number;
  [key: string]: unknown;
};

export type PathwayNetworkProps = {
  graph: PathwayGraph;
};

export const PathwayNetwork: React.FC<PathwayNetworkProps> = ({ graph }) => {
  const svgRef = useRef<SVGSVGElement | null>(null);
  const [showCofactors, setShowCofactors] = useState<boolean>(false);

  // Filter graph based on cofactor visibility
  const filteredGraph = useMemo<PathwayGraph>(() => {
    if (showCofactors) {
      return graph;
    }

    const keptNodes = graph.nodes.filter(
      n => !(n.type === "metabolite" && n.is_cofactor === true),
    );
    const keptIds = new Set(keptNodes.map(n => n.id));
    const keptLinks = graph.links.filter(
      l => keptIds.has(l.source) && keptIds.has(l.target),
    );

    return {
      ...graph,
      nodes: keptNodes,
      links: keptLinks,
    };
  }, [graph, showCofactors]);

  useEffect(() => {
    if (!svgRef.current) return;
    const svg = d3.select<SVGSVGElement, unknown>(svgRef.current);
    svg.selectAll("*").remove();

    const inner = svg.append("g");

    const zoom = d3
      .zoom<SVGSVGElement, unknown>()
      // upper bound stays constant, lower bound will be updated
      .scaleExtent([0.5, 7])
      .on("zoom", event => {
        inner.attr("transform", event.transform);
      });

    svg.call(zoom as any);

    const g = new dagreD3.graphlib.Graph()
      .setGraph({
        nodesep: 30,
        ranksep: 150,
        rankdir: "LR",
        marginx: 10,
        marginy: 10,
      })
      .setDefaultEdgeLabel(() => ({}));

    // Nodes
    filteredGraph.nodes.forEach(node => {
      if (node.type === "metabolite") {
        const container = document.createElement("div");
        container.id = `mol-svg-${node.id}`;
        container.style.width = `${chemWidth}px`;
        container.style.height = `${chemHeight}px`;

        g.setNode(node.id, {
          labelType: "html",
          label: container,
          width: chemWidth,
          height: chemHeight,
        });
      } else if (node.type === "reaction") {
        const labelText =
          node.name && node.name.trim().length > 0 ? node.name : node.id;
        g.setNode(node.id, {
          labelType: "string",
          label: labelText,
          width: chemWidth * 0.7,
          height: chemHeight * 0.4,
        });
      }
    });

    // Edges
    filteredGraph.links.forEach(link => {
      g.setEdge(link.source, link.target);
    });

    // Rounded corners
    g.nodes().forEach(v => {
      const n = g.node(v);
      n.rx = n.ry = 5;
    });

    const render = new (dagreD3 as any).render();
    render(inner as any, g as any);

    // Style overrides: white boxes, black borders, black arrows, Arial font
    inner
      .selectAll("g.node rect")
      .attr("fill", "#ffffff")
      .attr("stroke", "#000000")
      .attr("stroke-width", 1.2);

    inner
      .selectAll("g.node ellipse, g.node polygon")
      .attr("fill", "#ffffff")
      .attr("stroke", "#000000")
      .attr("stroke-width", 1.2);

    inner
      .selectAll("g.edgePath path")
      .attr("stroke", "#000000")
      .attr("stroke-width", 1.5)
      .attr("fill", "none");

    inner
      .selectAll("marker path")
      .attr("fill", "#000000");

    // Default node text styling
    inner
      .selectAll("g.node text")
      .attr("font-family", "Arial, sans-serif")
      .attr("font-size", 11);

    // Add metabolite labels
    filteredGraph.nodes.forEach(node => {
      if (node.type !== "metabolite") return;

      const labelText =
        node.name && node.name.trim().length > 0 ? node.name : node.id;

      const n = g.node(node.id);
      if (!n) return;

      // y position slightly above the top border so it does not sit on the line
      const topY = n.y - n.height / 3;
      const labelY = topY - 10;

      inner
        .append("text")
        .attr("x", n.x)
        .attr("y", labelY)
        .attr("text-anchor", "middle")
        .attr("font-size", 11)
        .attr("font-family", "Arial, sans-serif")
        .attr("fill", "#000")
        .attr("pointer-events", "none")
        .text(labelText);
    });

    // Fit and center, and set this scale as the minimum zoom
    const bounds = inner.node()?.getBBox();
    if (bounds && bounds.width && bounds.height) {
      const parent = svgRef.current.parentElement;
      if (!parent) return;

      const rect = parent.getBoundingClientRect();

      const midX = bounds.x + bounds.width / 2;
      const midY = bounds.y + bounds.height / 2;

      const initialScale =
        0.95 / Math.max(bounds.width / rect.width, bounds.height / rect.height);

      const translate: [number, number] = [
        rect.width / 2 - midX * initialScale,
        rect.height / 2 - midY * initialScale,
      ];

      // Use the initial scale as the minimum zoom so you cannot zoom out further
      zoom.scaleExtent([initialScale, 7]);

      const transform = d3.zoomIdentity
        .translate(translate[0], translate[1])
        .scale(initialScale);

      svg.transition().duration(300).call(zoom.transform as any, transform);
    }

    // Render SmilesSVG into metabolite nodes
    filteredGraph.nodes.forEach(node => {
      if (node.type !== "metabolite") return;
      const smiles = node.annotation?.smiles;
      if (!smiles) return;

      const svgContainer = document.getElementById(`mol-svg-${node.id}`);
      if (!svgContainer) return;

      const root = createRoot(svgContainer);
      root.render(
        <SmilesSVG smiles={smiles} width={chemWidth} height={chemHeight} />,
      );
    });
  }, [filteredGraph]);

  const handleDownloadSvg = () => {
    if (!svgRef.current) return;

    const svgNode = svgRef.current.cloneNode(true) as SVGSVGElement;
    svgNode.setAttribute("xmlns", "http://www.w3.org/2000/svg");
    svgNode.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");

    const serializer = new XMLSerializer();
    const source = serializer.serializeToString(svgNode);
    const blob = new Blob([source], {
      type: "image/svg+xml;charset=utf-8",
    });
    const url = URL.createObjectURL(blob);

    const a = document.createElement("a");
    a.href = url;
    a.download = `${graph.id ?? "pathway"}.svg`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  return (
    <div>
      <div
        style={{
          marginBottom: 8,
          display: "flex",
          alignItems: "center",
          gap: 16,
          flexWrap: "wrap",
        }}
      >
        <div>
          <span style={{ marginRight: 8 }}>Cofactors</span>
          <label style={{ marginRight: 12 }}>
            <input
              type="radio"
              name="cofactor-visibility"
              value="show"
              checked={showCofactors}
              onChange={() => setShowCofactors(true)}
              style={{ marginRight: 4 }}
            />
            Show
          </label>
          <label>
            <input
              type="radio"
              name="cofactor-visibility"
              value="hide"
              checked={!showCofactors}
              onChange={() => setShowCofactors(false)}
              style={{ marginRight: 4 }}
            />
            Hide
          </label>
        </div>

        <button
          type="button"
          onClick={handleDownloadSvg}
          style={{
            padding: "4px 10px",
            fontSize: 13,
            borderRadius: 4,
            border: "1px solid #ccc",
            backgroundColor: "#f7f7f7",
            cursor: "pointer",
          }}
        >
          Download SVG
        </button>
      </div>

      <div
        style={{
          width: "100%",
          height: "100%",
          borderRadius: 12,
          border: "1px solid #e0e0e0",
          padding: "12px 14px",
          boxShadow: "0 2px 4px rgba(0,0,0,0.04)",
          backgroundColor: "#ffffff",
          boxSizing: "border-box",
        }}
      >
        <svg
          ref={svgRef}
          style={{
            width: "100%",
            height: 500,
            display: "block",
          }}
        />
      </div>
    </div>
  );
};