# FastCollisionDetectionLib
C++ fast collision detection for uniform-distributed AABB particles using adaptive grid with implicit vectorization.

- 1000x speedup against naive brute-force algorithm for 40k particles, with uniform-distribution in range [0 - 1]).
- - 290x speedup when half of AABBs are 10x further than each other [0-1] and [10-11]
- - 230x speedup when half of AABBs are 10x far and a single AABB 100x far: [0-1] x N/2, [10-11] x N/2, [100-101] x1 
- Produced collision list does not contain duplicate pairs of collisions
- Particle data is not touched, work done only on pointers internally
- Currently it is adaptive, but needs optimizations on memory handling. 
- - On every cell-overflow, it stretches the cell to AABB of all particles and converts to a grid of 4x4x4 cells each with 4 capacity
- Implementation of IParticle is an AABB (axis-aligned bounding box) model 
- - In user defined particle (box as example here), methods (getMinX/Y/Z and getMaxX/Y/Z)  must return AABB corners of the underlying user-particle

```C++
#include"FastCollisionDetectionLib.h"

// any shape you like, but implements IParticle<float/double/etc> for the collision detection
class Box: public FastColDetLib::IParticle<float>
{
public:
	Box(float xPrm=0, float yPrm=0, float zPrm=0, float x2Prm=0, float y2Prm=0, float z2Prm=0, size_t idPrm=0)
	{
		x=xPrm;
		y=yPrm;
		z=zPrm;

		x2=x2Prm;
		y2=y2Prm;
		z2=z2Prm;
		id=idPrm;
	}

	const float getMaxX() const { return x>=x2?x:x2;}
	const float getMaxY() const { return y>=y2?y:y2;}
	const float getMaxZ() const { return z>=z2?z:z2;}

	const float getMinX() const { return x>=x2?x2:x;}
	const float getMinY() const { return y>=y2?y2:y;}
	const float getMinZ() const { return z>=z2?z2:z;}
	const int getId() const {return id;}

	~Box(){}

private:
	float x,y,z,x2,y2,z2;
	size_t id;
};

int main()
{
	const int n=40960;
	Box box[n];
	for(int i=0;i<n;i++)
	{
		auto x1 = /* generate your box coordinate */;
		auto y1 = /* generate your box coordinate */;
		auto z1 = /* generate your box coordinate */;

		// some computation to keep density same for different values of n
		float dx = 0.1f/std::pow(n,1.0f/3);
		box[i]=Box(x1,y1,z1,x1+dx,y1+dx,z1+dx,i);
	}

	FastColDetLib::AdaptiveGrid<float> grid(	64, // number of cells on X dimension
						64, // number of cells on Y dimension
						64, // number of cells on Z dimension
						35  // max particles per cell before cell becomes a sub-grid
	);
	FastColDetLib::BruteForce<float> bruteForce;
	std::vector<FastColDetLib::CollisionPair<float>> coll3D,coll3Dbrute;

	// adding pointers of all elements at once 
	// can be used multiple times, even from non-contiguous regions,
	// internal work uses pointers only (Iparticle<type>)
	grid.add(&box[0],n);
	bruteForce.add(&box[0],n);

	coll3D = grid.getCollisions();            // 14.5 milliseconds (fx8150 at 2.1GHz)
	coll3Dbrute = bruteForce.getCollisions(); // 15.3 seconds (fx8150 at 2.1GHz)

	for(auto& pair:coll3D)
	{
		int particleId1 = pair.getParticle1()->getId();
		int particleId2 = pair.getParticle2()->getId();
		// compute something using indices of particles that are colliding
	}

	return 0;
}
```
