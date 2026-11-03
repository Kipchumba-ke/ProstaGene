import { StrictMode, useEffect, useState } from "react";
import { createRoot } from "react-dom/client";
import "@fontsource/atkinson-hyperlegible/400.css";
import "@fontsource/atkinson-hyperlegible/700.css";
import "./styles.css";
import { useAppData } from "./data";
import { Header, Footer, MockBanner, type Route } from "./components/Layout";
import Home from "./pages/Home";
import Learn from "./pages/Learn";
import Clinician from "./pages/Clinician";
import Researcher from "./pages/Researcher";

const parseRoute = (): Route => {
  const h = window.location.hash.replace(/^#\/?/, "").split("?")[0];
  return (["learn", "clinician", "researcher"].includes(h) ? h : "home") as Route;
};

function App() {
  const [route, setRoute] = useState<Route>(parseRoute);
  const { data, error } = useAppData();
  useEffect(() => {
    const on = () => { setRoute(parseRoute()); window.scrollTo(0, 0); };
    window.addEventListener("hashchange", on);
    return () => window.removeEventListener("hashchange", on);
  }, []);
  useEffect(() => {
    const titles: Record<Route, string> = { home: "ProstaGene Explorer", learn: "Learn · ProstaGene", clinician: "Clinician · ProstaGene", researcher: "Researcher · ProstaGene" };
    document.title = titles[route];
  }, [route]);

  return (
    <>
      <a href="#main" className="sr-only" style={{ position: "absolute", left: -9999 }}>Skip to content</a>
      <Header route={route} />
      <main id="main">
        {error && <div className="wrap section"><div className="notice"><strong>Could not load data.</strong> {error}</div></div>}
        {!data && !error && <div className="wrap section"><p className="lead">Loading results…</p></div>}
        {data && (
          <>
            <MockBanner generated={data.summary.generated} />
            {route === "home" && <Home data={data} />}
            {route === "learn" && <Learn data={data} />}
            {route === "clinician" && <Clinician data={data} />}
            {route === "researcher" && <Researcher data={data} />}
          </>
        )}
      </main>
      <Footer generated={data?.summary.generated ?? ""} />
    </>
  );
}

createRoot(document.getElementById("root")!).render(<StrictMode><App /></StrictMode>);
