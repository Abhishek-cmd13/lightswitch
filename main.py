from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from scipy.spatial import Delaunay
import numpy as np

app = FastAPI()

class PolygonRequest(BaseModel):
    points: list[list[float]]

@app.post("/compute_centroids")
def compute_centroids(request: PolygonRequest):
    # Extract the points from the request body
    points = np.array(request.points)

    # Compute the Delaunay triangulation of the polygon
    tri = Delaunay(points)

    centroids = []

    # Compute the centroid of each triangle in the triangulation
    for t in tri.simplices:
        # Get the vertices of the current triangle
        tri_vertices = points[t]

        # Compute the centroid of the triangle
        centroid = np.mean(tri_vertices, axis=0)

        # Append the centroid to the list
        centroids.append(centroid.tolist())

    return {"centroids": centroids}

# Enable CORS (Cross-Origin Resource Sharing)
origins = [
    "http://localhost",
    "http://localhost:8000",
    "http://127.0.0.1:5500/index1.html"
    # Add more allowed origins as needed
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
