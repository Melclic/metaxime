// SmilesSVG.tsx
import { useEffect, useRef } from "react";
import SmilesDrawer from "smiles-drawer";

export type SmilesSVGProps = {
  smiles: string;
  width?: number;
  height?: number;
};

export function SmilesSVG({
  smiles,
  width = 300,
  height = 200
}: SmilesSVGProps) {
  //const containerRef = useRef<HTMLDivElement | null>(null);
  const containerRef = useRef<SVGSVGElement | null>(null);

  useEffect(() => {
    if (!smiles) return;
    if (!containerRef.current) return;

    // Clear previous SVG
    containerRef.current.innerHTML = "";

    const drawer = new SmilesDrawer.SvgDrawer({
      width,
      height,
      explicitHydrogens: false
    });

    SmilesDrawer.parse(
      smiles,
      (tree: any) => {
        // Draw directly into the container element
        drawer.draw(tree, containerRef.current, "light", false);
      },
      (error: any) => {
        console.error("Smiles parse error:", error);
      }
    );
  }, [smiles, width, height]);

  return (
    <svg
      ref={containerRef}
      data-smiles={smiles}
      width={width}
      height={height}
      style={{ display: "block" }}
    />
  );
}
