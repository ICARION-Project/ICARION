#!/usr/bin/env python3
"""ICARION Config Web GUI (local-only).

Run:
  pip install -r webui/requirements.txt
  python3 webui/app.py
Then open http://127.0.0.1:8000
"""

from pathlib import Path
from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi import Request

BASE_DIR = Path(__file__).resolve().parent
STATIC_DIR = BASE_DIR / "static"
REPO_ROOT = BASE_DIR.parent

app = FastAPI(title="ICARION Config Web GUI", docs_url=None, redoc_url=None)


@app.middleware("http")
async def disable_cache(request: Request, call_next):
  response = await call_next(request)
  response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
  response.headers["Pragma"] = "no-cache"
  response.headers["Expires"] = "0"
  return response

app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")


@app.get("/")
def index() -> FileResponse:
    return FileResponse(STATIC_DIR / "index.html")


@app.get("/api/load-json")
def load_json(path: str) -> JSONResponse:
  target = (REPO_ROOT / path).resolve()
  try:
    target.relative_to(REPO_ROOT)
  except ValueError as exc:
    raise HTTPException(status_code=400, detail="Path outside repo") from exc

  if not target.exists():
    raise HTTPException(status_code=404, detail="File not found")

  try:
    import json

    data = json.loads(target.read_text(encoding="utf-8"))
  except Exception as exc:
    raise HTTPException(status_code=400, detail=f"Invalid JSON: {exc}") from exc

  return JSONResponse(content=data)


if __name__ == "__main__":
    import uvicorn

    uvicorn.run("app:app", host="127.0.0.1", port=8000, reload=True)
