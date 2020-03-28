// Mini Game Engine
//
// Init() will fill the virtual world with 10.000.000 objects (represented by a point and an id 0:N)
// After Init() objects won't move anymore.
//
// Update() will update 10 observers, each represented as an AABB
//
// Render() will use those 10 observers AABBs to cull the existing objects, and output the count of objects overlapping the observers AABBs 


/*
		RESULTS:
			CPU: i7-6700HQ @ 2.60GHz
			Platform: x64

	|-----------------------|-----------------------|------------------------------------|
	|						| Initialization		|	Time Average time between frames |
	|	  Vanilla		    |	 3500.86ms			|	    1370.25ms					 |
	|-----------------------|-----------------------|------------------------------------|
	|		With			|						|									 |						
	|Octree implementation	|	 13712.6ms			|	     0.0386ms					 |
	------------------------|-----------------------|------------------------------------|

*/

#include "pch.h"
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <random>
#include <iostream>
#include <chrono>



///////////////////////////////////////////////////////////////////////////////
// Common types for fixed size over all platforms
typedef uint8_t uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;

typedef float  f32;
typedef double f64;

///////////////////////////////////////////////////////////////////////////////
// constants used throughout the program
 // enum : uint64 { OBJECT_COUNT   = 10000000};
enum : uint64 { OBJECT_COUNT = 10000000 };
enum : uint64 { OBSERVER_COUNT = 10 };
enum : uint32 { SEED = 12345 }; // hard-coded seed for deterministic random values



///////////////////////////////////////////////////////////////////////////////
// class to represent a 64 bit Vector
class Vect3
{
public:

	Vect3() { };
	explicit Vect3(f64 x) : m_x(x), m_y(x), m_z(x) {}
	Vect3(f64 x, f64 y, f64 z) : m_x(x), m_y(y), m_z(z) {}

	Vect3 operator+(const Vect3 &rOther) const { return Vect3(m_x + rOther.m_x, m_y + rOther.m_y, m_z + rOther.m_z); }
	Vect3 operator-(const Vect3 &rOther) const { return Vect3(m_x - rOther.m_x, m_y - rOther.m_y, m_z - rOther.m_z); }
	Vect3 operator*(f64 r) const { return Vect3(m_x*r, m_y*r, m_z*r); }
	Vect3& operator*=(float r) { m_x *= r; m_y *= r; m_z *= r;	return *this; }
	Vect3& operator+=(const Vect3& r) {	m_x += r.x(); m_y += r.y();	m_z += r.z(); return *this;	}

	f64 x() const { return m_x; }
	f64 y() const { return m_y; }
	f64 z() const { return m_z; }
	
	void setX(f64 x) { m_x = x; }
	void setY(f64 y) { m_y = y; }
	void setZ(f64 z) { m_z = z; }

private:
	f64 m_x, m_y, m_z;
};

///////////////////////////////////////////////////////////////////////////////
// class to represent an object in the virtual world
class CObject
{
public:
	CObject(){}
	CObject(const Vect3 &rPosition, uint64 nIdentifier) : m_vPosition(rPosition), m_nIdentifier(nIdentifier) {}

	const Vect3& GetPosition() const { return m_vPosition; }
	uint64 GetIdentifier() const { return m_nIdentifier; }

private:
	Vect3   m_vPosition;
	uint64 m_nIdentifier;
};

///////////////////////////////////////////////////////////////////////////////
// class to represent an Octree Node
class Octree {
	
	Vect3 m_Origin;         
	Vect3 m_HalfDimension;  //! Half the width/height/depth of this node

	
	Octree *m_Children[8]; 
	CObject *m_Object;   

	/* Child structure:

			child:	0 1 2 3 4 5 6 7
			x:      - - - - + + + +
			y:      - - + + - - + +
			z:      - + - + - + - +
	 */

public:

	Octree(const Vect3& origin, const Vect3& halfDimension)
		: m_Origin(origin), m_HalfDimension(halfDimension), m_Object(NULL) {
		
		for (int i = 0; i < 8; ++i)
			m_Children[i] = NULL;
	}
	
	~Octree() {
		
		for (int i = 0; i < 8; ++i)
			delete m_Children[i];
	}
	



	void insertObjectInTree(CObject* point) {

		if (isLeafNode()) {

			if (checkIfDoesntContainsObject()) {
				m_Object = point;
				return;
			}
			else {

				// Save Object that was here for a later re-insert
				CObject *oldPoint = m_Object;
				m_Object = NULL;

				splitNodeAndCreateNewOctTree();

				reinsertOldPointAndInsertNewPoint(oldPoint, point);

			}
		}
		else {

			//Recursive Insert Call
			insertInAppropiateChild(point);

		}
	}

	bool isLeafNode() const {

		return m_Children[0] == NULL;
	}

	bool checkIfDoesntContainsObject() {
		return m_Object == NULL;
	}

	void splitNodeAndCreateNewOctTree() {
		for (int i = 0; i < 8; ++i) {
			Vect3 newOrigin = m_Origin;
			newOrigin.setX(newOrigin.x() + m_HalfDimension.x() * (i & 4 ? .5f : -.5f));
			newOrigin.setY(newOrigin.y() + m_HalfDimension.y() * (i & 2 ? .5f : -.5f));
			newOrigin.setZ(newOrigin.z() + m_HalfDimension.z() * (i & 1 ? .5f : -.5f));
			m_Children[i] = new Octree(newOrigin, m_HalfDimension*.5f);
		}
	}

	void reinsertOldPointAndInsertNewPoint(CObject *oldPoint, CObject* point) {
		m_Children[getOctantContainingPoint(oldPoint->GetPosition())]->insertObjectInTree(oldPoint);
		m_Children[getOctantContainingPoint(point->GetPosition())]->insertObjectInTree(point);
	}

	void insertInAppropiateChild(CObject* point) {
		int octant = getOctantContainingPoint(point->GetPosition());
		m_Children[octant]->insertObjectInTree(point);
	}

	int getOctantContainingPoint(const Vect3& point) const {
		int oct = 0;
		if (point.x() >= m_Origin.x()) oct |= 4;
		if (point.y() >= m_Origin.y()) oct |= 2;
		if (point.z() >= m_Origin.z()) oct |= 1;
		return oct;
	}
	   	
	void getPointsContainedInBoundingBox(const Vect3& vMin, const Vect3& vMax, std::vector<uint64> & results) {
		
		if (isLeafNode()) {
			if (!checkIfDoesntContainsObject()) {
				const Vect3& p = m_Object->GetPosition();
				if (p.x() > vMax.x() || p.y() > vMax.y() || p.z() > vMax.z()) return;
				if (p.x() < vMin.x() || p.y() < vMin.y() || p.z() < vMin.z()) return;
				results.push_back(m_Object->GetIdentifier());
			}
		}
		else {

			for (int childOctree = 0; childOctree < 8; ++childOctree) {
				

				Vect3 maxCornerChild = m_Children[childOctree]->m_Origin + m_Children[childOctree]->m_HalfDimension;
				Vect3 minCornerChild = m_Children[childOctree]->m_Origin - m_Children[childOctree]->m_HalfDimension;

				// Check if the the Bounding Box is outside of the child's bounding box
				if (isOutOfBoundingBox(maxCornerChild,minCornerChild,vMin,vMax))
					continue;
				

				// this child is intersecting the box
				//recursive call
				m_Children[childOctree]->getPointsContainedInBoundingBox(vMin, vMax, results);
			}

		}

	}

	bool isOutOfBoundingBox(const Vect3 maxCornerChild, const Vect3 minCornerChild, const Vect3& vMin, const Vect3& vMax)
	{
		bool isOut = false;

		if (maxCornerChild.x() < vMin.x() || maxCornerChild.y() < vMin.y() || maxCornerChild.z() < vMin.z()) isOut = true;
		if (minCornerChild.x() > vMax.x() || minCornerChild.y() > vMax.y() || minCornerChild.z() > vMax.z()) isOut = true;

		return isOut;
	}

};


///////////////////////////////////////////////////////////////////////////////
// class to represent an Axis-Aligned-Bounding-Box (AABB)
class AABB
{
public:
	AABB() : m_vMin(0), m_vMax(0) {}
	AABB(const Vect3 &rMin, const Vect3 &rMax) : m_vMin(rMin) , m_vMax(rMax) {}

	bool Overlaps(const Vect3 &rPoint) const
	{
		if(rPoint.x() < m_vMin.x()) return false;
		if(rPoint.y() < m_vMin.y()) return false;
		if(rPoint.z() < m_vMin.z()) return false;
		if(rPoint.x() > m_vMax.x()) return false;
		if(rPoint.y() > m_vMax.y()) return false;
		if(rPoint.z() > m_vMax.z()) return false;

		return true;
	}
	Vect3 m_vMin;
	Vect3 m_vMax;
private:
	
};

///////////////////////////////////////////////////////////////////////////////
// class to represent an observer into the virtual world
class CObserver
{
public:
	CObserver() : m_Bounds() {}

	void UpdateAABB(const AABB &rBounds) { m_Bounds = rBounds; }

	AABB GetBounds() const { return m_Bounds; }
private:
	AABB m_Bounds;
};

///////////////////////////////////////////////////////////////////////////////
struct WorldState
{
	//A Mersenne Twister pseudo-random generator of 64-bit numbers with a state size of 19937 bits.
	std::mt19937_64         rnd_gen;

	std::vector<CObject*>   arrObjects;
	std::vector<CObserver*> arrObservers;
	Octree *octree;
};


///////////////////////////////////////////////////////////////////////////////
void Init(WorldState &rWorldState, uint32 nSeed)
{
	auto start = std::chrono::high_resolution_clock::now();

	rWorldState.rnd_gen.seed(nSeed);
	std::uniform_real_distribution<f64> dist(-1e10, +1e10);	
	
	// Create a new Octree in the origin	
	rWorldState.octree = new Octree(Vect3(0, 0, 0), Vect3(1e10, 1e10, 1e10));


	std::cout << "Building Octree" << std::endl;
	for (uint64 i = 0; i < OBJECT_COUNT; ++i)
	{
				
		const Vect3 vPosition = Vect3(dist(rWorldState.rnd_gen), dist(rWorldState.rnd_gen), dist(rWorldState.rnd_gen));

		CObject* newObject = new CObject(vPosition, i);

		rWorldState.arrObjects.push_back(newObject);		
		rWorldState.octree->insertObjectInTree(newObject);

	}
		
	std::cout << "Inserted points to Octree" << std::endl;
	
	
	for (uint64 i = 0; i < OBSERVER_COUNT; ++i)
	{
		rWorldState.arrObservers.push_back(new CObserver());
	}

	std::cout << "Init took " << std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - start).count() << " ms" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
f64 Update(WorldState &rWorldState)
{
	auto start = std::chrono::high_resolution_clock::now();

	std::uniform_real_distribution<f64> dist_center(-1e10, +1e10);	
	std::uniform_real_distribution<f64> dist_extends(+20000000, +40000000);	

	for (CObserver *pObserver : rWorldState.arrObservers)
	{
		const Vect3 vCenter = Vect3(dist_center(rWorldState.rnd_gen), dist_center(rWorldState.rnd_gen), dist_center(rWorldState.rnd_gen));
		const f64 fExtends = dist_extends(rWorldState.rnd_gen);

		pObserver->UpdateAABB(AABB(vCenter - Vect3(fExtends), vCenter + Vect3(fExtends)));
	}

	return std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - start).count();
}

///////////////////////////////////////////////////////////////////////////////
void checkOverlappingObjectsWithObserver(CObserver* observer,std::vector<uint64> &rOverlappingObjects, WorldState &rWorldState) {
	auto start = std::chrono::high_resolution_clock::now();

	
	std::vector<CObject*> results;
	rWorldState.octree->getPointsContainedInBoundingBox(observer->GetBounds().m_vMin, observer->GetBounds().m_vMax, rOverlappingObjects);
	
	auto T = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - start).count();
	
}

///////////////////////////////////////////////////////////////////////////////
f64 Render(WorldState &rWorldState, std::vector<uint64> &rOverlappingObjects)
{
	auto start = std::chrono::high_resolution_clock::now();

	
	for (CObserver *pObserver : rWorldState.arrObservers)
	{
		
		checkOverlappingObjectsWithObserver(pObserver,rOverlappingObjects, rWorldState);
	}

	//RAW
	/*for (CObserver *pObserver : rWorldState.arrObservers)
	{
		for (CObject* pObject : rWorldState.arrObjects)
		{
			if (pObserver->GetBounds().Overlaps(pObject->GetPosition()))
			{
				rOverlappingObjects.push_back(pObject->GetIdentifier());
			}
		}
	}*/
	
	return std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - start).count();
}

///////////////////////////////////////////////////////////////////////////////
int main()
{
	WorldState state;
	int oct = 0;
	

	Init(state, SEED);
	
	
	for (size_t i = 0; i < 10; ++i)
	{
		std::vector<uint64> arrOverlappingObjects;

		const f64 fTimeUpdate = Update(state);
		const f64 fTimeRender = Render(state, arrOverlappingObjects);

		std::cout << "Frame " << i << ". Update: " << fTimeUpdate << " ms. Render: " << fTimeRender << " ms. Visible Objects: " << arrOverlappingObjects.size() << " [";
		
		for(uint64 nObject : arrOverlappingObjects)
			std::cout << nObject << " ";

		std::cout << "]" << std::endl;
	}

	std::cin.ignore();


	//// object tear-down/shutdown code can be omitted (in fact the OS will just clean that up for us more efficiently)
 	return EXIT_SUCCESS;
}

