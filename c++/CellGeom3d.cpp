#include "CellGeom3d.h"

// Superclass constructor is also called below.
CellGeom3d::CellGeom3d(int idNew, char geomTypeNew, 
	std::vector<int> vertIndicesNew, std::vector<Vertex*> vertPointersNew)
	: CellGeomBase(idNew, geomTypeNew, vertIndicesNew, vertPointersNew)
{
	assignFaceVertIndices();
	computeFaceCentrePoints();
	computeFaceAreas();
	computeVolume();
}

void CellGeom3d::setNeighIndices(std::vector<int> neighIndicesNew)
{
	neighIndices = neighIndicesNew;
}

void CellGeom3d::setNeighPointers(std::vector<CellGeom3d*> neighPointersNew)
{
	neighPointers = neighPointersNew;
	computeDistancesBetweenCells();
}

void CellGeom3d::assignFaceVertIndices()
{
	// Hexahedre.
	if (geomType == 12) {
		std::vector<int> faceTop{ 4,5,6,7 };
		faceVertIndices.push_back(faceTop);
		std::vector<int> faceBack{ 2,3,7,6 };
		faceVertIndices.push_back(faceBack);
		std::vector<int> faceLeft{ 3,0,4,7 };
		faceVertIndices.push_back(faceLeft);
		std::vector<int> faceFront{ 0,1,5,4 };
		faceVertIndices.push_back(faceFront);
		std::vector<int> faceRight{ 1,2,6,5 };
		faceVertIndices.push_back(faceRight);
		std::vector<int> faceBottom{ 3,2,1,0 };
		faceVertIndices.push_back(faceBottom);
	}
	// Wedge.
	else if (geomType == 13) {
		std::vector<int> faceTop{ 3,4,5 };
		faceVertIndices.push_back(faceTop);
		std::vector<int> face_side_01{ 0,3,5,2 };
		faceVertIndices.push_back(face_side_01);
		std::vector<int> face_side_02{ 1,4,3,0 };
		faceVertIndices.push_back(face_side_02);
		std::vector<int> face_side_03{ 2,5,4,0 };
		faceVertIndices.push_back(face_side_03);
		std::vector<int> faceBottom{ 0,2,1 };
		faceVertIndices.push_back(faceBottom);
	}
}

void CellGeom3d::computeFaceCentrePoints()
{
	for (size_t i = 0; i < faceVertIndices.size(); i++) {
		Vertex faceCentrePoint;
		faceCentrePoint.x = 0.0;
		faceCentrePoint.y = 0.0;
		faceCentrePoint.z = 0.0;
		for (size_t j = 0; j < faceVertIndices.at(i).size(); j++) {
			faceCentrePoint.x += vertPointers.at(faceVertIndices.at(i).at(j))->x;
			faceCentrePoint.y += vertPointers.at(faceVertIndices.at(i).at(j))->y;
			faceCentrePoint.z += vertPointers.at(faceVertIndices.at(i).at(j))->z;
		}
		if (faceVertIndices.at(i).size() > 0) {
			faceCentrePoint.x /= faceVertIndices.at(i).size();
			faceCentrePoint.y /= faceVertIndices.at(i).size();
			faceCentrePoint.z /= faceVertIndices.at(i).size();
		}
		faceCentrePoints.push_back(faceCentrePoint);
	}
}

void CellGeom3d::computeDistancesBetweenCells()
{
	for (size_t i = 0; i < faceCentrePoints.size(); i++) {
		if (neighPointers.at(i) != 0){
			Vertex centrePointNeigh = neighPointers.at(i)->getCentrePoint();
			Vertex vector0 = createVector(centrePoint, faceCentrePoints.at(i));
			double distance0 = computeVectorLength(vector0);
			Vertex vector1 = createVector(faceCentrePoints.at(i), 
				centrePointNeigh);
			double distance1 = computeVectorLength(vector1);
			distances.push_back(distance0 + distance1);
		}
		else {
			distances.push_back(0.0);
		}
	}
}

void CellGeom3d::computeFaceAreas()
{
	// Hexahedre.
	if (geomType == 12) {
		for (size_t i = 0; i < faceVertIndices.size(); i++) {
			double faceArea;
			Vertex vector1 = createVector(*vertPointers.at(faceVertIndices.at(i).at(0)),
				*vertPointers.at(faceVertIndices.at(i).at(3)));
			Vertex vector2 = createVector(*vertPointers.at(faceVertIndices.at(i).at(0)),
				*vertPointers.at(faceVertIndices.at(i).at(1)));
			Vertex crossProduct = computeVectorCrossProduct(vector1, vector2);
			faceArea = 0.5 * computeVectorLength(crossProduct);
			vector1 = createVector(*vertPointers.at(faceVertIndices.at(i).at(2)),
				*vertPointers.at(faceVertIndices.at(i).at(1)));
			vector2 = createVector(*vertPointers.at(faceVertIndices.at(i).at(2)),
				*vertPointers.at(faceVertIndices.at(i).at(3)));
			crossProduct = computeVectorCrossProduct(vector1, vector2);
			faceArea += 0.5 * computeVectorLength(crossProduct);
			faceAreas.push_back(faceArea);
		}
	}
	// Wedge.
	else if (geomType == 13) {
		// Top side.
		Vertex vector1 = createVector(*vertPointers.at(faceVertIndices.at(0).at(0)),
			*vertPointers.at(faceVertIndices.at(0).at(2)));
		Vertex vector2 = createVector(*vertPointers.at(faceVertIndices.at(0).at(0)),
			*vertPointers.at(faceVertIndices.at(0).at(1)));
		Vertex crossProduct = computeVectorCrossProduct(vector1, vector2);
		faceAreas.push_back(0.5 * computeVectorLength(crossProduct));

		// Side faces.
		for (size_t i = 1; i < faceVertIndices.size() - 1; i++) {
			double faceArea;
			Vertex vector1 = createVector(*vertPointers.at(faceVertIndices.at(i).at(0)),
				*vertPointers.at(faceVertIndices.at(i).at(3)));
			Vertex vector2 = createVector(*vertPointers.at(faceVertIndices.at(i).at(0)),
				*vertPointers.at(faceVertIndices.at(i).at(1)));
			Vertex crossProduct = computeVectorCrossProduct(vector1, vector2);
			faceArea = 0.5 * computeVectorLength(crossProduct);
			vector1 = createVector(*vertPointers.at(faceVertIndices.at(i).at(2)),
				*vertPointers.at(faceVertIndices.at(i).at(1)));
			vector2 = createVector(*vertPointers.at(faceVertIndices.at(i).at(2)),
				*vertPointers.at(faceVertIndices.at(i).at(3)));
			crossProduct = computeVectorCrossProduct(vector1, vector2);
			faceArea += 0.5 * computeVectorLength(crossProduct);
			faceAreas.push_back(faceArea);
		}

		// Bottom side.
		vector1 = createVector(*vertPointers.at(faceVertIndices.at(4).at(0)),
			*vertPointers.at(faceVertIndices.at(4).at(2)));
		vector2 = createVector(*vertPointers.at(faceVertIndices.at(4).at(0)),
			*vertPointers.at(faceVertIndices.at(4).at(1)));
		crossProduct = computeVectorCrossProduct(vector1, vector2);
		faceAreas.push_back(0.5 * computeVectorLength(crossProduct));
	}
}

void CellGeom3d::computeVolume()
{
	// Hexahedre.
	if (geomType == 12) {
		Vertex vector = createVector(faceCentrePoints.at(5), 
			faceCentrePoints.at(0));
		volume = computeVectorLength(vector) * 0.5 * 
			(faceAreas.at(5) + faceAreas.at(0));
	}
	// Wedge.
	else if (geomType == 13) {
		Vertex vector = createVector(faceCentrePoints.at(4),
			faceCentrePoints.at(0));
		volume = computeVectorLength(vector) * 0.5 *
			(faceAreas.at(4) + faceAreas.at(0));
	}
}
