# Elizabeth

\## Project goal

Estimate the Elizabeth line’s impact on London retail rateable values. Outcome: Δlog(RV)=ln(RV2023)−ln(RV2017). Treatment: distance/buffers (0–400 m, 400–800 m) to nearest Elizabeth line station. Controls: baseline PTAL (2015), IMD 2019, CAZ/Town Centre flags, and Borough fixed effects.



\## Folder structure

\- data/raw/        Raw inputs (do not modify)

\- data/interim/    Intermediate artefacts (reprojection, filters)

\- data/clean/      Analysis-ready datasets

\- docs/            Documentation (data inventory, methods)

\- logs/            Cleaning log (`cleaning\_log.md`)

\- outputs/         Figures and tables



\## Naming \& CRS conventions

\- File naming: `agency\_topic\_year\[\_extra].ext` (e.g., `tfl\_elizabeth\_line\_stations\_2022.geojson`)

\- Primary CRS: \*\*EPSG:27700 (British National Grid)\*\*; WGS84 (EPSG:4326) for web maps.

\- Core fields: `property\_id, rv\_2017, rv\_2023, use\_class, borough, ptal, imd\_score, caz\_flag, town\_centre\_flag, dist\_elizabeth, ring\_0\_400, ring\_400\_800`.



\## Data sources

See `docs/data\_inventory.md` and `docs/data\_inventory.csv` for links, licenses, and usage notes. Keep original archives in `data/raw/`. Write any derived files to `data/interim/` or `data/clean/` and record changes in `logs/cleaning\_log.md`.



