import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// Relative base: the built site works from any sub-path (e.g. /<repo>/explorer/ on GitHub Pages).
// Routing uses the URL hash, so no server rewrites are needed.
export default defineConfig({
  plugins: [react()],
  base: process.env.EXPLORER_BASE ?? "./",
  build: {
    outDir: process.env.EXPLORER_OUT ?? "dist",
    emptyOutDir: true,
  },
});
