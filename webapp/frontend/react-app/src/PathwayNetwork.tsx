import React, { useEffect, useRef } from "react";
import { createRoot } from "react-dom/client";
// @ts-ignore
import * as d3 from "d3";
import * as dagreD3 from 'dagre-d3-es';
import { SmilesSVG } from "./SmilesSVG";

type MetaboliteNode = {
  id: string;
  type: "metabolite";
  name: string;
  topology: string;
  annotation?: {
    smiles?: string;
    [key: string]: unknown;
  };
  [key: string]: unknown;
};

type ReactionNode = {
  id: string;
  type: "reaction";
  name: string;
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
  width?: number | string;
  height?: number | string;
};

/**
 * Draws a metabolic network using dagre-d3, similar to the old makeNetwork / drawNetwork
 * implementation, but using the SmilesSVG React component to render metabolite structures.
 *
 * - Metabolites are rendered as SVG structures (via SmilesSVG inside HTML labels / foreignObject)
 * - Reactions are rendered as simple rounded boxes with their id as text
 */
export const PathwayNetwork: React.FC<PathwayNetworkProps> = ({
  graph,
  width = "100%",
  height = 500,
}) => {
  const svgRef = useRef<SVGSVGElement | null>(null);

  useEffect(() => {
    if (!svgRef.current) return;
    const svg = d3.select<SVGSVGElement, unknown>(svgRef.current);
    svg.selectAll("*").remove();

    const inner = svg.append("g");

    const zoom = d3
      .zoom<SVGSVGElement, unknown>()
      .scaleExtent([0.5, 7])
      .on("zoom", (event: d3.D3ZoomEvent<SVGSVGElement, unknown>) => {
        inner.attr("transform", event.transform);
      });

    svg.call(zoom as any);

    // Create dagre graph
    const g = new dagreD3.graphlib.Graph()
      .setGraph({
        nodesep: 30,
        ranksep: 150,
        rankdir: "LR",
        marginx: 10,
        marginy: 10,
      })
      .setDefaultEdgeLabel(() => ({}));

    const chemWidth = 140;
    const chemHeight = 100;

    // Nodes
    graph.nodes.forEach(node => {
      if (node.type === "metabolite") {
        // HTML container that will host the SmilesSVG React component
        const container = document.createElement("div");
        container.id = `mol-${node.id}`;
        container.style.width = `${chemWidth}px`;
        container.style.height = `${chemHeight}px`;

        g.setNode(node.id, {
          labelType: "html",
          label: container,
          width: chemWidth,
          height: chemHeight,
        });
      } else if (node.type === "reaction") {
        g.setNode(node.id, {
          labelType: "string",
          label: node.id,
          width: chemWidth * 0.6,
          height: chemHeight * 0.4,
        });
      }
    });

    // Edges
    graph.links.forEach(link => {
      g.setEdge(link.source, link.target);
    });

    // Rounded corners
    g.nodes().forEach(v => {
      const n = g.node(v);
      n.rx = n.ry = 5;
    });

    // Render
    const render = new (dagreD3 as any).render();
    render(inner as any, g as any);

    // Fit / center the graph (similar to zoomFit)
    const bounds = inner.node()?.getBBox();
    if (bounds && bounds.width && bounds.height) {
      const parent = svgRef.current.parentElement;
      const fullWidth = parent?.clientWidth ?? 800;
      const fullHeight = parent?.clientHeight ?? 400;

      const midX = bounds.x + bounds.width / 2;
      const midY = bounds.y + bounds.height / 2;

      const scale =
        0.85 / Math.max(bounds.width / fullWidth, bounds.height / fullHeight);
      const translate: [number, number] = [
        fullWidth / 2 - midX * scale,
        fullHeight / 2 - midY * scale,
      ];

      const transform = d3.zoomIdentity
        .translate(translate[0], translate[1])
        .scale(scale);

      svg.transition().duration(300).call(zoom.transform as any, transform);
    }

    // Render SmilesSVG inside each metabolite container
    graph.nodes.forEach(node => {
      if (node.type !== "metabolite") return;
      const smiles = node.annotation?.smiles;
      if (!smiles) return;

      const container = document.getElementById(`mol-${node.id}`);
      if (!container) return;

      const root = createRoot(container);
      root.render(
        <SmilesSVG smiles={smiles} width={chemWidth} height={chemHeight} />
      );
    });
  }, [graph, height, width]);

  return <svg ref={svgRef} width={width} height={height} />;
};