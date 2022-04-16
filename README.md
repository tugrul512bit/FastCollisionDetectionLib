# FastCollisionDetectionLib
C++ fast collision detection for uniform(and non-uniform)-distributed AABB particles using adaptive grid with implicit vectorization. 

## Non-Sparse Adaptive Grid

- 10 million dynamic particles AABB collision check per second against static grid of 8000 particles
- 1000x speedup against naive brute-force algorithm for 40k particles (static vs static), with uniform-distribution in range [0 - 1]).
- - teapot-in-stadium problem is partially solved by "adaptive" grid:
- - 290x speedup when half of AABBs are 10x further than each other [0-1] and [10-11]
- - 230x speedup when half of AABBs are 10x far and a single AABB 100x far: [0-1] x N/2, [10-11] x N/2, [100-101] x1 
- Produced collision list does not contain duplicate pairs of collisions
- Particle data is not touched, work done only on pointers internally
- Currently it is adaptive, but needs optimizations on memory handling. 
- - On every cell-overflow, it stretches the cell to AABB of all particles and converts to a grid of 4x4x4 cells each with 4 capacity
- Implementation of IParticle is an AABB (axis-aligned bounding box) model 
- - In user defined particle (box as example here), methods (getMinX/Y/Z and getMaxX/Y/Z)  must return AABB corners of the underlying user-particle

## Sparse - Linear - Adaptive Grid

- More than 2000x speedup against naive brute-force, in single-thread for 10000 particles
- 60 FPS for 20000 particles and ~25000 collision pairs on 2.1GHz FX8150 single-thread
- - 100+ FPS for new CPUs
- Better performance stability compared to non-sparse version
- Better SIMD support on all-pairs computation method using tiled-computing
- Non-zero based object-id values supported (getId() method in IParticle<float> interface)

```C++
// prepare memory pool
FastColDetLib::MemoryPool memPool;

// map grid to a volume of cube between corners of (0,0,0) and (10005,10005,10005)
FastColDetLib::AdaptiveGridV2 grid(memPool,0,0,0,10005,10005,10005);

// implement IParticle<float>
struct AABBofPointCloud: public FastColDetLib::IParticle<float>
{
   ...
   	const CoordType getMaxX()const {return xmax;}
	const CoordType getMaxY()const {return ymax;}
	const CoordType getMaxZ()const {return zmax;}
	const CoordType getMinX()const {return xmin;}
	const CoordType getMinY()const {return ymin;}
	const CoordType getMinZ()const {return zmin;}
	const int getId()const {return id;}
    ...
};

// initialize AABB vector
std::vector<AABBofPointCloud> AABBs;

while(simulation)
{
    // clear tree data
    grid.clear();
    
    // add particles that implement IParticle<float> into grid
    grid.addParticles(N,AABBs.data());
    
    // build tree
    grid.buildTree();
    
    // compute all-pairs collision array 
    // 60FPS on FX8150 2.1GHz single-thread for 20000 particles with less than 29000 collisions
    // 100FPS on Xeon Gold 5215 2.9GHz single-thread for 20000 particles with less than 29000 collisions
    std::vector<std::pair<int,int>> vec = grid2_0.findCollisionsAll();
    
    // the vec contains id-values of particles that have their AABBs collide so that you can do further fine-grained collision checks between them
}

```

## Multithreaded Tree of Sparse - Linear - Adaptive Grid

- 60 FPS for 40000 particles' AABB all-pair computations
- Bottlenecked by RAM bandwidth and mutex-array locking throughput
- Only zero-based object-id values supported
- Work load is balanced on particles, not volumes, this makes better distribution of particles on threads but causes duplicated work due to merging of results from all leaf nodes

```C++
    // 7 threads with load-balancing by a tree, mapped to (0,0,0)-(10005,10005,10005) region
    AdaptiveGridTree<7> test(0,0,0,10005,10005,10005);

    for(int i=0;i<100;i++)
    {
		size_t nano;
		{
			FastColDetLib::Bench bench(&nano);
			
			// clear contents of memory pool of adaptive grid tree
			test.clear();
			
			// AABBs is a vector of objects that implements IParticle<float> interface 
	                // (only zero-based object-id values supported from getId() method of IParticle<float>)
			test.addParticles(N,AABBs.data()); 
			
			// non-duplicate pairs of collisions 
			const auto coll = test.computeAllPairs();
			std::cout<<"c="<<coll.size()<<std::endl;
		}
		std::cout<<"t="<<nano<<std::endl;
    }
```



For details, please visit ![https://github.com/tugrul512bit/FastCollisionDetectionLib/wiki](https://github.com/tugrul512bit/FastCollisionDetectionLib/wiki) wiki page.

Working demo for non-sparse (old version) adaptive grid (requires linking pthread for header and gomp/fopenmp for this demo):

```C++
#include"FastCollisionDetectionLib.h"

// user-object, with IParticle<coordinate_type> interface for querying AABB within the API
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

#include<iostream>
#include<omp.h>
#include"Generator.h"
int main()
{

	oofrng::Generator<64> gen;

	{
		// d^3 number of particles
		const int d = 20;
		const int n=d*d*d;
		std::cout<<"n="<<n<<std::endl;
		std::vector<Box> box(n);
		for(int i=0;i<n;i++)
		{
			auto x1 = gen.generate1Float()*d;
			auto y1 = gen.generate1Float()*d;
			auto z1 = gen.generate1Float()*d;

			float dx = 0.05f;
			float dx2 = 0.05f;
			box[i]=Box(x1,y1,z1,x1+dx2+dx*gen.generate1Float(),y1+dx2+dx*gen.generate1Float(),z1+dx2+dx*gen.generate1Float(),i);
		}

		// thread-pool releases itself once it is out of scope
		// single instance of thread pool can be used for multiple adaptive grids
		FastColDetLib::ThreadPool<float> thr;
		FastColDetLib::AdaptiveGrid<float> grid(thr,-1,-1,-1,d+1,d+1,d+1);
		FastColDetLib::BruteForce<float> bruteForce;
		std::vector<FastColDetLib::CollisionPair<float>> coll3D,coll3Dbrute;


		std::cout<<"add"<<std::endl;
		bruteForce.add(&box[0],n);
		std::cout<<"compute"<<std::endl;
		bool bfor = true;


		for(int k=0;k<40;k++)
		{
			size_t nano1,nano2;

			{
				{
					FastColDetLib::Bench bench(&nano1);
					{
						FastColDetLib::Bench bench(&nano1);
						grid.clear();
					}
					std::cout<<"grid-static-object-clear: "<<nano1/1000000000.0<<" s "<<std::endl;
					{
						FastColDetLib::Bench bench(&nano1);
						grid.add(box.data(),n);
					}
					std::cout<<"grid-static-object-add ("<<n<<" particles AABB): "<<nano1/1000000000.0<<" s "<<std::endl;
					{
						FastColDetLib::Bench bench(&nano1);
						coll3D = grid.getCollisions();
					}
					std::cout<<"grid-static-object-compute: "<<nano1/1000000000.0<<" s "<<coll3D.size()<<std::endl;
				}
				std::cout<<"grid-static-total: "<<nano1/1000000000.0<<" s "<<coll3D.size()<<std::endl;
				std::vector<float> rand(1000000*3);
				gen.generate(rand.data(),1000000*3);
				{
					std::mutex mut;
					std::vector<FastColDetLib::IParticle<float>*> res;
					{
						FastColDetLib::Bench bench(&nano1);
						#pragma omp parallel for
						for(int j=0;j<100;j++)
						{
							std::vector<FastColDetLib::IParticle<float>*> resTmp;
							for(int i=0;i<1000;i++)
							{
								auto x = rand[i*3]*d;
								auto y = rand[i*3+1]*d;
								auto z = rand[i*3+2]*d;
								auto item = Box(x,y,z,x+0.25,y+0.25,z+0.25);
								auto collisions = grid.getDynamicCollisionListFor(&item);
								std::copy(collisions.begin(),collisions.end(),std::back_inserter(resTmp));
							}
							std::lock_guard<std::mutex> lg(mut);
							std::copy(resTmp.begin(),resTmp.end(),std::back_inserter(res));
						}

					}
					std::cout<<"grid-compute-dynamic (100k particles AABB): "<<nano1/1000000000.0<<" s "<<res.size()<<std::endl;
				}

			}


			if(bfor)
			{
				FastColDetLib::Bench bench(&nano2);
				coll3Dbrute = bruteForce.getCollisions();
			}
			if(bfor)
				std::cout<<"Brute-force ("<<n<<" particles AABB): "<<nano2/1000000000.0<<" s "<<coll3Dbrute.size()<<std::endl;
			std::cout<<"------------------------------------------------------------------------------------"<<std::endl;
		}


		if(bfor)
		{
			if(coll3D.size() != coll3Dbrute.size())
			{
				std::cout<<"ERROR: not equal sizes of collision lists"<<std::endl;
				std::cout<<coll3D.size()<<" vs "<< coll3Dbrute.size()<<std::endl;
				const int sz = std::min(coll3D.size(),coll3Dbrute.size());
				const int sz2 = sz<10?sz:10;
				for(int i=0;i<sz2;i++)
				{
					if((coll3D[i].getParticle1()->getId()!=coll3Dbrute[i].getParticle1()->getId()) &&
							(coll3D[i].getParticle2()->getId()!=coll3Dbrute[i].getParticle2()->getId())
							)
					{
						std::cout<<"ERRRROOOR!"<<std::endl;
						std::cout<<coll3D[i].getParticle1()->getId()<<"<-->"<<coll3D[i].getParticle2()->getId()<<"        "<<coll3Dbrute[i].getParticle1()->getId()<<"<-->"<<coll3Dbrute[i].getParticle2()->getId()<<std::endl;
					}
				}
			}
			else
			{
				const size_t sz = coll3D.size();
				for(size_t i=0;i<sz;i++)
				{

					if((coll3D[i].getParticle1()->getId()!=coll3Dbrute[i].getParticle1()->getId()) &&
							(coll3D[i].getParticle2()->getId()!=coll3Dbrute[i].getParticle2()->getId())
							)
					{
						std::cout<<"ERRRROOOR!"<<std::endl;
						std::cout<<coll3D[i].getParticle1()->getId()<<"<-->"<<coll3D[i].getParticle2()->getId()<<"        "<<coll3Dbrute[i].getParticle1()->getId()<<"<-->"<<coll3Dbrute[i].getParticle2()->getId()<<std::endl;
					}
				}
			}
		}
		std::cout<<"--"<<std::endl;


	}
	return 0;
}

```

output for FX8150 at 2.1GHz:

```
n=8000
add
compute
grid-static-object-clear: 6.4614e-05 s 
grid-static-object-add (8000 particles AABB): 0.00656692 s 
grid-static-object-compute: 0.00195991 s 16
grid-static-total: 0.00866529 s 16
grid-compute-dynamic (100k particles AABB): 0.0137853 s 3600
Brute-force (8000 particles AABB): 0.541069 s 16
------------------------------------------------------------------------------------
grid-static-object-clear: 0.00143487 s 
grid-static-object-add (8000 particles AABB): 0.00523593 s 
grid-static-object-compute: 0.00182263 s 16
grid-static-total: 0.00853913 s 16
grid-compute-dynamic (100k particles AABB): 0.00835323 s 4300
Brute-force (8000 particles AABB): 0.541662 s 16
------------------------------------------------------------------------------------
grid-static-object-clear: 0.00161331 s 
grid-static-object-add (8000 particles AABB): 0.00521717 s 
grid-static-object-compute: 0.00180544 s 16
grid-static-total: 0.0086723 s 16
grid-compute-dynamic (100k particles AABB): 0.00864333 s 3200
Brute-force (8000 particles AABB): 0.54128 s 16
------------------------------------------------------------------------------------
grid-static-object-clear: 0.00151348 s 
grid-static-object-add (8000 particles AABB): 0.00524639 s 
grid-static-object-compute: 0.00180144 s 16
grid-static-total: 0.00860804 s 16
grid-compute-dynamic (100k particles AABB): 0.010513 s 3200
Brute-force (8000 particles AABB): 0.541246 s 16
------------------------------------------------------------------------------------
grid-static-object-clear: 0.00160501 s 
grid-static-object-add (8000 particles AABB): 0.00521408 s 
grid-static-object-compute: 0.00176137 s 16
grid-static-total: 0.00862669 s 16
grid-compute-dynamic (100k particles AABB): 0.0120599 s 3000
Brute-force (8000 particles AABB): 0.54391 s 16
------------------------------------------------------------------------------------
grid-static-object-clear: 0.00174257 s 
grid-static-object-add (8000 particles AABB): 0.00548767 s 
grid-static-object-compute: 0.00182023 s 16
grid-static-total: 0.00909577 s 16
grid-compute-dynamic (100k particles AABB): 0.0083235 s 3500
Brute-force (8000 particles AABB): 0.54521 s 16
------------------------------------------------------------------------------------

```

This is ~8500x performance for checking a particle collision against a grid of static objects and  ~50x performance for static-static collision check (or dynamic-dynamic).

More particles:

```
n=64000
add
compute
grid-static-object-clear: 1.3058e-05 s 
grid-static-object-add (64000 particles AABB): 0.0347275 s 
grid-static-object-compute: 0.0106954 s 114
grid-static-total: 0.0455527 s 114
grid-compute-dynamic (100k particles AABB): 0.0154652 s 3700
Brute-force (64000 particles AABB): 35.473 s 114

```

777x performance for static vs static collision checking, more than 2000x performance for static vs dynamic.
