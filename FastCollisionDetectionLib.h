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
#include<queue>
#include<thread>
#include<mutex>
#include<set>
#include<functional>
#include<condition_variable>

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


	template<typename CoordType>
	class AdaptiveGrid;


	using GridDataType = int;

	// Fixed grid of cells (adaptive if a cell overflows)
	template<typename CoordType>
	class FixedGridFields
	{
	public:
		FixedGridFields(const int w, const int h, const int d,	const int s,
				const CoordType minXp, const CoordType minYp, const CoordType minZp,
				const CoordType maxXp, const CoordType maxYp, const CoordType maxZp):
					width(w),height(h),depth(d),storage(s),stride(storage+1 /* w*h*d*/),
					minX(minXp),minY(minYp),minZ(minZp),maxX(maxXp),maxY(maxYp),maxZ(maxZp)
		{

			grid.resize(width*height*depth*(storage+1));
		}

		~FixedGridFields()
		{

		}

		inline
		void clearGrid() noexcept
		{
			for(int k=0;k<depth;k++)
			{
				const int ki = k*width*height;
				for(int j=0;j<height;j++)
				{
					const int jki = j*width + ki;
					for(int i=0;i<width;i++)
					{
						//grid[i+jki]=0; // 0: no particle, -1: a grid
						grid[(i+jki)*stride]=0; // 0: no particle, -1: a grid
					}
				}
			}
		}

		inline
		const GridDataType numParticlesAtCell(const int indexX, const int indexY, const int indexZ) const noexcept
		{
			//return grid[indexX + indexY*width + indexZ*width*height];
			return grid[(indexX + indexY*width + indexZ*width*height)*stride];
		}



		inline
		const GridDataType numParticlesAtCell(const int indexX, const int indexY, const int indexZ, int& computedCellIndex) const noexcept
		{
			//computedCellIndex=indexX + indexY*width + indexZ*width*height;
			computedCellIndex=(indexX + indexY*width + indexZ*width*height)*stride;
			return grid[computedCellIndex];
		}

		inline
		void numParticlesAtCellIncrement(const int computedCellIndex) noexcept
		{
			grid[computedCellIndex]++;
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
		void setCellData(const int cellId, const int laneId, const GridDataType data) noexcept
		{
			//grid[cellId + (laneId+1)*stride] = data; // +1: first item is number of particles or indicator of another grid (adaptiveness)
			grid[cellId + (laneId+1)] = data; // +1: first item is number of particles or indicator of another grid (adaptiveness)
		}

		inline
		const GridDataType getCellData(const int cellId, const int laneId) const noexcept
		{
			//return grid[cellId + (laneId+1)*stride]; // +1: first item is number of particles or indicator of another grid (adaptiveness)
			return grid[cellId + (laneId+1)]; // +1: first item is number of particles or indicator of another grid (adaptiveness)
		}

		std::vector<IParticle<CoordType>*> particles;

		std::map<IParticle<CoordType>*,std::map<IParticle<CoordType>*,bool>> coll;

		// first integer in cell = number of particles. -1: another grid (adaptiveness)
		std::vector<GridDataType> grid;
		std::map<IParticle<CoordType>*,std::map<IParticle<CoordType>*,bool>> mapping;
		std::vector<CollisionPair<CoordType>> result;
		const int width;
		const int height;
		const int depth;
		const int storage;
		const int stride;
		const CoordType minX;
		const CoordType minY;
		const CoordType minZ;
		const CoordType maxX;
		const CoordType maxY;
		const CoordType maxZ;

	};


	template<typename CoordType>
	struct Cmd
	{
		AdaptiveGrid<CoordType>* grid;
		std::mutex* mut;
		std::map<IParticle<CoordType>*,std::map<IParticle<CoordType>*,bool>>* mapping;
		bool* completed;
	};



	template<typename T>
	class SyncQueue
	{
	public:
		SyncQueue(){}
		void push(T t)
		{
			std::unique_lock<std::mutex> lc(m);
			q.push(t);
			c.notify_one();
		}

		void push2(T t)
		{
			std::unique_lock<std::mutex> lc(m);
			q.push(t);
			c.notify_all();
		}

		T pop()
		{
			std::unique_lock<std::mutex> lc(m);
			while(q.empty())
			{
				c.wait(lc);
			}
			T result = q.front();
			q.pop();
			return result;
		}
	private:
		std::queue<T> q;
		std::mutex m;
		std::condition_variable c;
	};

	class ThreadPoolDestructor
	{
	public:
		ThreadPoolDestructor(std::function<void(void)>  dd)
		{
			destroyer = dd;
		}
		std::function<void(void)> destroyer;
		~ThreadPoolDestructor()
		{
			destroyer();
		}
	};

	template<typename CoordType>
	class ThreadPool
	{
	public:
		ThreadPool():ctr(std::make_shared<int>())
		{
			worker=std::make_shared<std::vector<std::thread>>();
			mut=std::make_shared<std::mutex>();
			msg=std::make_shared<std::vector<int>>();
			q=std::make_shared<std::vector<std::shared_ptr<SyncQueue<Cmd<CoordType>>>>>();
			destro=std::make_shared<ThreadPoolDestructor>([&,this](){this->stop();});
			for(int i=0;i<7;i++)
			{
				q->push_back(std::make_shared<SyncQueue<Cmd<CoordType>>>());
				msg->push_back(1);
				worker->push_back(std::thread(
						[&,i]()
						{
							bool work = true;
							while(work)
							{
								{
									{
										std::lock_guard<std::mutex> lg(*mut);
										work=((*msg)[i]>0);
									}

									Cmd<CoordType> cmd = (*q)[i]->pop();
									if(cmd.grid==nullptr)
										break;
									auto collisions = cmd.grid->getCollisions();

									{
										std::lock_guard<std::mutex> lg(*cmd.mut);
										for(auto& c:collisions)
										{
											(*cmd.mapping)[c.getParticle1()][c.getParticle2()]=true;
										}
										*cmd.completed=true;
									}

								}
							}
						}
				));
			}
		}


		void compute(Cmd<CoordType> cmd)
		{
			(*q)[(*ctr)++%7]->push(cmd);
		}

		void stop()
		{
			for(unsigned int i=0;i<worker->size();i++)
			{
				std::lock_guard<std::mutex> lg(*mut);
				(*msg)[i]=0;
				Cmd<CoordType> cmd;
				cmd.grid=nullptr;
				(*q)[i]->push2(cmd);
			}

			for(unsigned int i=0;i<worker->size();i++)
			{

				(*worker)[i].join();
			}
		}
	private:
		std::shared_ptr<int> ctr;
		std::shared_ptr<std::vector<std::thread>> worker;
		std::shared_ptr<std::mutex> mut;
		std::shared_ptr<std::vector<int>> msg;
		std::shared_ptr<std::vector<std::shared_ptr<SyncQueue<Cmd<CoordType>>>>> q;
		std::shared_ptr<ThreadPoolDestructor> destro;
	};

	template<typename CoordType>
	class AdaptiveGrid
	{

protected:

		AdaptiveGrid(ThreadPool<CoordType> thr, int depthPrm,
				const CoordType minX, const CoordType minY, const CoordType minZ,
				const CoordType maxX, const CoordType maxY, const CoordType maxZ):thrPool(thr)
		{
			depth = std::make_shared<int>();

			*depth=depthPrm;

			fields = std::make_shared<FixedGridFields<CoordType>>(2,2,2,8+*depth*2,minX,minY,minZ,maxX,maxY,maxZ);
			subGrid = std::make_shared<std::vector<AdaptiveGrid<CoordType>>>();
			fields->clearGrid();
		}

		// compute collision between given particle and the already-prepared static object grid (after add(..) and getCollisions(..))
		std::set<IParticle<CoordType>*> getDynamicCollisionSetFor(IParticle<CoordType>* particle)
		{
			std::set<IParticle<CoordType>*> result;
			// AABB box of particle
			const CoordType minx = particle->getMinX();
			const CoordType miny = particle->getMinY();
			const CoordType minz = particle->getMinZ();

			const CoordType maxx = particle->getMaxX();
			const CoordType maxy = particle->getMaxY();
			const CoordType maxz = particle->getMaxZ();

			const CoordType xDim = fields->maxX - fields->minX;
			const CoordType yDim = fields->maxY - fields->minY;
			const CoordType zDim = fields->maxZ - fields->minZ;

			const int w = fields->getWidth();
			const int h = fields->getHeight();
			const int d = fields->getDepth();
			const CoordType stepX = xDim/w;
			const CoordType stepY = yDim/h;
			const CoordType stepZ = zDim/d;


			const int cellIndexX = std::floor((minx - fields->minX)/stepX);
			const int cellIndexY = std::floor((miny - fields->minY)/stepY);
			const int cellIndexZ = std::floor((minz - fields->minZ)/stepZ);

			const int cellIndexX2 = std::floor((maxx - fields->minX)/stepX);
			const int cellIndexY2 = std::floor((maxy - fields->minY)/stepY);
			const int cellIndexZ2 = std::floor((maxz - fields->minZ)/stepZ);

			for(int zz = cellIndexZ; zz<=cellIndexZ2; zz++)
				for(int yy = cellIndexY; yy<=cellIndexY2; yy++)
					for(int xx = cellIndexX; xx<=cellIndexX2; xx++)
					{
						if(xx<0 || yy<0 || zz<0 || xx>=w || yy>=h || zz>=d)
							continue;
						// num particles in cell
						int computedCellIndex = 0;

						const int n = fields->numParticlesAtCell(xx, yy, zz, computedCellIndex);

						// if this is a cell
						if(n>=0)
						{
							for(int j=0;j<n;j++)
							{
								const int k = fields->getCellData(computedCellIndex,j);



								const CoordType minx2 = fields->particles[k]->getMinX();
								const CoordType maxx2 = fields->particles[k]->getMaxX();
								if(intersectDim(minx, maxx, minx2, maxx2))
								{

									const CoordType miny2 = fields->particles[k]->getMinY();
									const CoordType maxy2 = fields->particles[k]->getMaxY();
									if(intersectDim(miny, maxy, miny2, maxy2))
									{

										const CoordType minz2 = fields->particles[k]->getMinZ();
										const CoordType maxz2 = fields->particles[k]->getMaxZ();
										if(intersectDim(minz, maxz, minz2, maxz2))
										{
											result.insert(fields->particles[k]);
										}
									}
								}


							}

						}
						else // if this is a grid
						{

							const int gridIdx = -n-1;


							auto subResult = (*subGrid)[gridIdx].getDynamicCollisionListFor(particle);
							result.insert(subResult.begin(),subResult.end());
						}
					}

			return result;
		}

public:
		AdaptiveGrid(ThreadPool<CoordType> thr,
				const CoordType minX, const CoordType minY, const CoordType minZ,
				const CoordType maxX, const CoordType maxY, const CoordType maxZ):thrPool(thr)
		{
			depth = std::make_shared<int>();

			*depth=0;

			fields = std::make_shared<FixedGridFields<CoordType>>(16,16,16,2,minX,minY,minZ,maxX,maxY,maxZ);
			subGrid = std::make_shared<std::vector<AdaptiveGrid<CoordType>>>();
			fields->clearGrid();
		}





		void clear()
		{
			fields->clearGrid();
			subGrid->clear();
			fields->particles.clear();
		}

		// add static particle object pointers to compute all-vs-all comparisons in an optimized way
		// the generated internal data is also used for static vs dynamic collision checking
		template<typename Derived>
		void add(Derived * particlesPrm, int numParticlesToAdd)
		{

			const int w = fields->getWidth();
			const int h = fields->getHeight();
			const int d = fields->getDepth();

			const int sto = fields->getStorage();
			// grid
			const CoordType xDim = fields->maxX - fields->minX;
			const CoordType yDim = fields->maxY - fields->minY;
			const CoordType zDim = fields->maxZ - fields->minZ;

			// cell
			const CoordType stepX = xDim/w;
			const CoordType stepY = yDim/h;
			const CoordType stepZ = zDim/d;

			for(int ii=0;ii<numParticlesToAdd;ii++)
			{
				const int i=fields->particles.size();

				// AABB box of particle
				const CoordType minx = (particlesPrm+ii)->getMinX();
				const CoordType miny = (particlesPrm+ii)->getMinY();
				const CoordType minz = (particlesPrm+ii)->getMinZ();

				const CoordType maxx = (particlesPrm+ii)->getMaxX();
				const CoordType maxy = (particlesPrm+ii)->getMaxY();
				const CoordType maxz = (particlesPrm+ii)->getMaxZ();





				const int cellIndexX = std::floor((minx - fields->minX)/stepX);
				const int cellIndexY = std::floor((miny - fields->minY)/stepY);
				const int cellIndexZ = std::floor((minz - fields->minZ)/stepZ);

				const int cellIndexX2 = std::floor((maxx - fields->minX)/stepX);
				const int cellIndexY2 = std::floor((maxy - fields->minY)/stepY);
				const int cellIndexZ2 = std::floor((maxz - fields->minZ)/stepZ);
				bool toCell = false;
				// "gather" operations on neighbor cells should be cache-friendly
				for(int zz = cellIndexZ; zz<=cellIndexZ2; zz++)
					for(int yy = cellIndexY; yy<=cellIndexY2; yy++)
						for(int xx = cellIndexX; xx<=cellIndexX2; xx++)
						{
							if(xx<0 || yy<0 || zz<0 || xx>=w || yy>=h || zz>=d)
								continue;
							// num particles in cell
							int computedCellIndex = 0;

							const int n = fields->numParticlesAtCell(xx, yy, zz, computedCellIndex);

							// if this is a cell
							if(n>=0)
							{
								// if capacity not reached
								if(n<sto)
								{

									fields->setCellData(computedCellIndex,n,i);
									fields->numParticlesAtCellIncrement(computedCellIndex);
									toCell=true;
								}
								else
								{


									const int idx = subGrid->size();

									AdaptiveGrid<CoordType> newGrid(thrPool,*depth+1,fields->minX+stepX*xx,fields->minY+stepY*yy,fields->minZ+stepZ*zz,
											fields->minX+(stepX)*(xx+1),fields->minY+(stepY)*(yy+1),fields->minZ+(stepZ)*(zz+1));
									subGrid->push_back(newGrid);

									fields->setCellData(computedCellIndex,-1,-(idx+1));



									for(int ins = 0; ins<n; ins++)
									{
										newGrid.add(fields->particles[fields->getCellData(computedCellIndex,ins)],1);
									}



									newGrid.add(particlesPrm+ii,1);
								}
							}
							else // if this is a grid
							{

								const int gridIdx = -n-1;


								(*subGrid)[gridIdx].add(particlesPrm+ii,1);

							}
						}

				if(toCell)
					fields->particles.push_back(particlesPrm+ii);

			}


		}




		inline
		const bool intersectDim(const CoordType minx, const CoordType maxx, const CoordType minx2, const CoordType maxx2) const noexcept
		{
			return !((maxx < minx2) || (maxx2 < minx));
		}

		// compute collision between given particle and the already-prepared static object grid (after add(..) and getCollisions(..))
		std::vector<IParticle<CoordType>*> getDynamicCollisionListFor(IParticle<CoordType>* particle)
		{
			std::set<IParticle<CoordType>*> result;
			// AABB box of particle
			const CoordType minx = particle->getMinX();
			const CoordType miny = particle->getMinY();
			const CoordType minz = particle->getMinZ();

			const CoordType maxx = particle->getMaxX();
			const CoordType maxy = particle->getMaxY();
			const CoordType maxz = particle->getMaxZ();

			const CoordType xDim = fields->maxX - fields->minX;
			const CoordType yDim = fields->maxY - fields->minY;
			const CoordType zDim = fields->maxZ - fields->minZ;

			const int w = fields->getWidth();
			const int h = fields->getHeight();
			const int d = fields->getDepth();
			const CoordType stepX = xDim/w;
			const CoordType stepY = yDim/h;
			const CoordType stepZ = zDim/d;


			const int cellIndexX = std::floor((minx - fields->minX)/stepX);
			const int cellIndexY = std::floor((miny - fields->minY)/stepY);
			const int cellIndexZ = std::floor((minz - fields->minZ)/stepZ);

			const int cellIndexX2 = std::floor((maxx - fields->minX)/stepX);
			const int cellIndexY2 = std::floor((maxy - fields->minY)/stepY);
			const int cellIndexZ2 = std::floor((maxz - fields->minZ)/stepZ);

			for(int zz = cellIndexZ; zz<=cellIndexZ2; zz++)
				for(int yy = cellIndexY; yy<=cellIndexY2; yy++)
					for(int xx = cellIndexX; xx<=cellIndexX2; xx++)
					{
						if(xx<0 || yy<0 || zz<0 || xx>=w || yy>=h || zz>=d)
							continue;
						// num particles in cell
						int computedCellIndex = 0;

						const int n = fields->numParticlesAtCell(xx, yy, zz, computedCellIndex);

						// if this is a cell
						if(n>=0)
						{
							for(int j=0;j<n;j++)
							{
								const int k = fields->getCellData(computedCellIndex,j);



								const CoordType minx2 = fields->particles[k]->getMinX();
								const CoordType maxx2 = fields->particles[k]->getMaxX();
								if(intersectDim(minx, maxx, minx2, maxx2))
								{

									const CoordType miny2 = fields->particles[k]->getMinY();
									const CoordType maxy2 = fields->particles[k]->getMaxY();
									if(intersectDim(miny, maxy, miny2, maxy2))
									{

										const CoordType minz2 = fields->particles[k]->getMinZ();
										const CoordType maxz2 = fields->particles[k]->getMaxZ();
										if(intersectDim(minz, maxz, minz2, maxz2))
										{
											result.insert(fields->particles[k]);
										}
									}
								}


							}

						}
						else // if this is a grid
						{

							const int gridIdx = -n-1;


							auto subResult = (*subGrid)[gridIdx].getDynamicCollisionSetFor(particle);
							result.insert(subResult.begin(),subResult.end());
						}
					}
			std::vector<IParticle<CoordType>*> resultVec;
			resultVec.insert(resultVec.end(),result.begin(),result.end());
			return resultVec;
		}

		// returns collision pairs between static objects (and prepares internal data for future dynamic object collision checking), ordered
		std::vector<CollisionPair<CoordType>> getCollisions()
		{
			fields->mapping.clear();
			fields->result.clear();


			const int w = fields->getWidth();
			const int h = fields->getHeight();
			const int d = fields->getDepth();


			// check neighbor cells for a collision of another AABB particle


			std::mutex mut;
			bool completed[8]={false,false,false,false,false,false,false,false};
			int completedCtr = 0;
			// "gather" operations on neighbor cells should be cache-friendly
			for(int zz = 0; zz<d; zz++)
				for(int yy = 0; yy<h; yy++)
					for(int xx = 0; xx<w; xx++)
			{


				int computedCellIndex0 = 0;
				const int n = fields->numParticlesAtCell(xx, yy, zz, computedCellIndex0);

				// if a cell
				if(n>=0)
				{
					const int computedCellIndex = computedCellIndex0;
					if(n<2)
						continue;

					for(int j=0;j<n-1;j++)
					{
						const int k = fields->getCellData(computedCellIndex,j);
						for(int i=j+1;i<n;i++)
						{

							const int k2 = fields->getCellData(computedCellIndex,i);
							if(fields->particles[k]->getId()<fields->particles[k2]->getId())
							{
								const CoordType minx = fields->particles[k]->getMinX();
								const CoordType maxx = fields->particles[k]->getMaxX();
								const CoordType minx2 = fields->particles[k2]->getMinX();
								const CoordType maxx2 = fields->particles[k2]->getMaxX();
								if(intersectDim(minx, maxx, minx2, maxx2))
								{
									const CoordType miny = fields->particles[k]->getMinY();
									const CoordType maxy = fields->particles[k]->getMaxY();
									const CoordType miny2 = fields->particles[k2]->getMinY();
									const CoordType maxy2 = fields->particles[k2]->getMaxY();
									if(intersectDim(miny, maxy, miny2, maxy2))
									{
										const CoordType minz = fields->particles[k]->getMinZ();
										const CoordType maxz = fields->particles[k]->getMaxZ();
										const CoordType minz2 = fields->particles[k2]->getMinZ();
										const CoordType maxz2 = fields->particles[k2]->getMaxZ();
										if(intersectDim(minz, maxz, minz2, maxz2))
										{
											std::lock_guard<std::mutex> lg(mut);
											fields->mapping[fields->particles[k]][fields->particles[k2]]=true;
										}
									}
								}
							}
						}
					}
				}
				else // if a grid
				{
					// if at specific layer, enable threads
					if(*depth==8)
					{
						Cmd<CoordType> cmd;
						cmd.completed=&completed[completedCtr++];
						cmd.mapping=&fields->mapping;
						cmd.mut=&mut;
						cmd.grid=subGrid->data()-n-1;
						thrPool.compute(cmd);
					}
					else
					{
						auto collisions = (subGrid->data()-n-1)->getCollisions();

						{

							for(auto& c:collisions)
							{
								fields->mapping[c.getParticle1()][c.getParticle2()]=true;
							}

						}
					}
				}
			}

			bool comp = false;
			while(!comp)
			{
				comp=true;
				{
					std::lock_guard<std::mutex> lg(mut);
					for(int cmd=0;cmd<completedCtr;cmd++)
					{
						comp&=completed[cmd];
					}
				}
				std::this_thread::yield();
			}
			for(auto& m:fields->mapping)
			{
				for(auto& m2:m.second)
				{

					fields->result.push_back(CollisionPair<CoordType>(m.first,m2.first));
				}
			}




			return fields->result;
		}
private:
		std::shared_ptr<FixedGridFields<CoordType>> fields;
		std::shared_ptr<std::vector<AdaptiveGrid<CoordType>>> subGrid;
		ThreadPool<CoordType> thrPool;
		std::shared_ptr<int> depth;

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
