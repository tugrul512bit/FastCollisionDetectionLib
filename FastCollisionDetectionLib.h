/*
 * FastCollisionDetectionLib.h
 *
 *  Created on: Mar 11, 2022
 *      Author: tugrul
 */

#ifndef FASTCOLLISIONDETECTIONLIB_H_
#define FASTCOLLISIONDETECTIONLIB_H_

#include<algorithm>
#include<vector>
#include<map>
#include<chrono>
#include<memory>
#include<math.h>
#include<iostream>

namespace FastColDetLib
{
	class Bench
	{
	public:
		Bench(size_t * targetPtr)
		{
			target=targetPtr;
			t1 =  std::chrono::duration_cast< std::chrono::nanoseconds >(std::chrono::high_resolution_clock::now().time_since_epoch());
		}

		~Bench()
		{
			t2 =  std::chrono::duration_cast< std::chrono::nanoseconds >(std::chrono::high_resolution_clock::now().time_since_epoch());
			*target= t2.count() - t1.count();
		}
	private:
		size_t * target;
		std::chrono::nanoseconds t1,t2;
	};

	/*
	 * interface to build various objects that can collide each other
	 *
	 */

	template<typename CoordType>
	class IParticle
	{
	public:
		virtual const CoordType getMaxX()const =0;
		virtual const CoordType getMaxY()const =0;
		virtual const CoordType getMaxZ()const =0;

		virtual const CoordType getMinX()const =0;
		virtual const CoordType getMinY()const =0;
		virtual const CoordType getMinZ()const =0;
		virtual const int getId()const =0;

		const bool intersectX(IParticle<CoordType>* p)
		{
			return !((getMaxX() < p->getMinX()) || (p->getMaxX() < getMinX()));
		}

		const bool intersectY(IParticle<CoordType>* p)
		{
			return !((getMaxY() < p->getMinY()) || (p->getMaxY() < getMinY()));
		}

		const bool intersectZ(IParticle<CoordType>* p)
		{
			return !((getMaxZ() < p->getMinZ()) || (p->getMaxZ() < getMinZ()));
		}

		virtual ~IParticle(){};
	};

	template<typename CoordType>
	class CollisionPair
	{
	public:
		CollisionPair(IParticle<CoordType>* p1Prm=nullptr, IParticle<CoordType>* p2Prm=nullptr)
		{
			p1=p1Prm;
			p2=p2Prm;
		}

		IParticle<CoordType>* getParticle1() const
		{
			return p1;
		}

		IParticle<CoordType>* getParticle2() const
		{
			return p2;
		}
	private:
		IParticle<CoordType> * p1;
		IParticle<CoordType> * p2;
	};




	// Fixed grid of cells (adaptive if a cell overflows)
	template<typename CoordType>
	class FixedGridFields
	{
	public:
		FixedGridFields(const int w=64, const int h=64, const int d=64, const int s=63):width(w),height(h),depth(d),storage(s),stride(w*h*d)
		{
			grid.resize(4096 + width*height*depth*(storage+1));
		}

		inline
		void clearGrid() noexcept
		{
			for(int k=0;k<depth;k++)
				for(int j=0;j<height;j++)
					for(int i=0;i<width;i++)
					{
						grid[i+j*width+k*width*height]=0; // 0: no particle, -1: a grid
					}
		}

		inline
		const int numParticlesAtCell(const int indexX, const int indexY, const int indexZ) const noexcept
		{
			return grid[indexX + indexY*width + indexZ*width*height];
		}



		inline
		const int numParticlesAtCell(const int indexX, const int indexY, const int indexZ, int& computedCellIndex) const noexcept
		{
			computedCellIndex=indexX + indexY*width + indexZ*width*height;
			return grid[computedCellIndex];
		}

		inline
		const int numParticlesAtCellIncrement(const int indexX, const int indexY, const int indexZ, int& computedCellIndex) noexcept
		{
			computedCellIndex=indexX + indexY*width + indexZ*width*height;
			return grid[computedCellIndex]++;
		}

		inline
		const int getWidth () const noexcept { return width;};

		inline
		const int getHeight () const noexcept { return height;};

		inline
		const int getDepth () const noexcept { return depth;};

		inline
		const int getStorage () const noexcept { return storage;};

		inline
		void setCellData(const int cellId, const int laneId, const int data) noexcept
		{
			grid[cellId + (laneId+1)*stride] = data; // +1: first item is number of particles or indicator of another grid (adaptiveness)
		}

		inline
		const int getCellData(const int cellId, const int laneId) const noexcept
		{
			return grid[cellId + (laneId+1)*stride]; // +1: first item is number of particles or indicator of another grid (adaptiveness)
		}


	private:
		// first integer in cell = number of particles. -1: another grid (adaptiveness)
		std::vector<int> grid;
		std::vector<char> scanIndexList;
		const int width;
		const int height;
		const int depth;
		const int storage;
		const int stride;
	};

	template<typename CoordType>
	class FixedGrid
	{
public:
		// 63 storage = 64 int per cell (1 for count)
		FixedGrid(const int width=64, const int height=64, const int depth=64, const int storage=63)
		{
			fields = std::make_shared<FixedGridFields<CoordType>>(width,height,depth,storage);
		}

		template<typename Derived>
		void add(Derived * particlesPrm, int numParticlesToAdd)
		{
			for(int i=0;i<numParticlesToAdd;i++)
			{
				particles.push_back(static_cast<IParticle<CoordType>*>(particlesPrm+i));
			}
		}

		std::vector<CollisionPair<CoordType>> getCollisions()
		{
			std::vector<CollisionPair<CoordType>> result;
			CoordType minX = particles[0]->getMinX();
			CoordType minY = particles[0]->getMinY();
			CoordType minZ = particles[0]->getMinZ();
			CoordType maxX = particles[0]->getMaxX();
			CoordType maxY = particles[0]->getMaxY();
			CoordType maxZ = particles[0]->getMaxZ();
			const int sz = particles.size();
			for(int i=0;i<sz;i++)
			{
				const IParticle<CoordType>* const p = particles[i];
				const CoordType minx = p->getMinX();
				const CoordType maxx = p->getMaxX();

				const CoordType miny = p->getMinY();
				const CoordType maxy = p->getMaxY();

				const CoordType minz = p->getMinZ();
				const CoordType maxz = p->getMaxZ();

				if(minx<minX)
					minX = minx;

				if(miny<minY)
					minY = miny;

				if(minz<minZ)
					minZ = minz;

				if(maxx>maxX)
					maxX = maxx;

				if(maxy>maxY)
					maxY = maxy;

				if(maxz>maxZ)
					maxZ = maxz;
			}

			CoordType deltax = maxX - minX;
			CoordType deltay = maxY - minY;
			CoordType deltaz = maxZ - minZ;

			// 2% widen to allow centers of particles on edges
			minX -= deltax*0.01f;
			minY -= deltay*0.01f;
			minZ -= deltaz*0.01f;
			maxX += deltax*0.01f;
			maxY += deltay*0.01f;
			maxZ += deltaz*0.01f;

			// new delta
			const int w = fields->getWidth();
			const int h = fields->getHeight();
			const int d = fields->getDepth();
			const CoordType deltaX = (maxX - minX)/w;
			const CoordType deltaY = (maxY - minY)/h;
			const CoordType deltaZ = (maxZ - minZ)/d;



			const int maxParticlesPerCell = fields->getStorage();

			fields->clearGrid();
			for(int i=0;i<sz;i++)
			{
				const IParticle<CoordType>* const p = particles[i];
				const CoordType minx = p->getMinX();
				const CoordType maxx = p->getMaxX();

				const CoordType miny = p->getMinY();
				const CoordType maxy = p->getMaxY();

				const CoordType minz = p->getMinZ();
				const CoordType maxz = p->getMaxZ();



				const CoordType stepX = (minx - minX);
				const CoordType stepY = (miny - minY);
				const CoordType stepZ = (minz - minZ);

				const CoordType stepX2 = (maxx - minX);
				const CoordType stepY2 = (maxy - minY);
				const CoordType stepZ2 = (maxz - minZ);

				const int idxX = std::floor(stepX/deltaX);
				const int idxY = std::floor(stepY/deltaY);
				const int idxZ = std::floor(stepZ/deltaZ);
				const int idxX2 = std::floor(stepX2/deltaX);
				const int idxY2 = std::floor(stepY2/deltaY);
				const int idxZ2 = std::floor(stepZ2/deltaZ);

				const int indexX = idxX % w;
				const int indexY = idxY % h;
				const int indexZ = idxZ % d;
				const int indexX2 = idxX2 % w;
				const int indexY2 = idxY2 % h;
				const int indexZ2 = idxZ2 % d;

				for(int kk=indexZ;kk<=indexZ2;kk++)
					for(int jj=indexY;jj<=indexY2;jj++)
						for(int ii=indexX;ii<=indexX2;ii++)
						{
							int cellIndex=-1;
							const int indexParticleInCell = fields->numParticlesAtCellIncrement(ii,jj,kk,cellIndex);

							if(indexParticleInCell<maxParticlesPerCell)
							{
								// add to all intersecting cells for fast query
								fields->setCellData(cellIndex,indexParticleInCell,i);
							}
							else
							{
								std::cout<<"Error: Grid-cell overflow. Max particle: "<<maxParticlesPerCell<<" added: "<<(indexParticleInCell+1)<<"  cell-id:"<<cellIndex<<std::endl;
							}
						}

			}

			// cell-based computation (SIMD-vectorized) in steps of 16
			std::map<int,std::map<int,bool>> coll;
			for(int k=0;k<d;k++)
				for(int j=0;j<h;j++)
					for(int i=0;i<w;i++)
					{
						int cellIndex=-1;
						const int n=fields->numParticlesAtCell(i,j,k,cellIndex);

						// brute-force within cell, n is small so its better than brute force for N

						if(n<2)
							continue;
						for(int pt = 0; pt<n-1; pt++)
						{
							const int idx = fields->getCellData(cellIndex,pt);
							for(int pt2 = pt+1; pt2<n; pt2++)
							{
								const int idx2 = fields->getCellData(cellIndex,pt2);

								if(		particles[idx]->intersectX(particles[idx2]) &&
										particles[idx]->intersectY(particles[idx2]) &&
										particles[idx]->intersectZ(particles[idx2])
										)
								{
									coll[idx][idx2]=true;
								}
							}


						}
					}

			for(auto& c:coll)
			{
				for(auto& c2:c.second)
				{
					result.push_back(CollisionPair<CoordType>(particles[c.first],particles[c2.first]));
				}
			}
			return result;
		}
private:
		std::vector<IParticle<CoordType>*> particles;
		std::shared_ptr<FixedGridFields<CoordType>> fields;
	};


	// axis-aligned bounding-box collision detection
	template<typename CoordType>
	class BruteForce
	{
public:
		BruteForce()
		{

		}

		template<typename Derived>
		void add(Derived * particlesPrm, int numParticlesToAdd)
		{
			for(int i=0;i<numParticlesToAdd;i++)
				particles.push_back(static_cast<IParticle<CoordType>*>(particlesPrm+i));
		}

		std::vector<CollisionPair<CoordType>> getCollisions()
		{
			std::vector<CollisionPair<CoordType>> result;
			idMap.clear();
			collisionPairs.clear();
			const int sz = particles.size();
			for(int i=0;i<sz-1;i++)
			{
				idMap[particles[i]->getId()]=particles[i];

				for(int j=i+1;j<sz;j++)
				{

					if( particles[i]->intersectX(particles[j]) && particles[i]->intersectY(particles[j]) && particles[i]->intersectZ(particles[j]))
					{
						collisionPairs.push_back(CollisionPair<CoordType>(particles[i],particles[j]));
					}
				}
			}
			std::sort(collisionPairs.begin(),collisionPairs.end(),[](CollisionPair<CoordType>& c1, CollisionPair<CoordType>& c2){
				return c1.getParticle1()->getId()<c2.getParticle1()->getId();
			});
			result=collisionPairs;
			return result;
		}
private:
		std::vector<IParticle<CoordType>*> particles;
		std::vector<CollisionPair<CoordType>> collisionPairs;
		std::map<int,IParticle<CoordType>*> idMap;
	};
}



#endif /* FASTCOLLISIONDETECTIONLIB_H_ */
