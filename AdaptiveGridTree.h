/*
 * AdaptiveGridTree.h
 *
 *  Created on: Apr 14, 2022
 *      Author: tugrul
 */

#ifndef ADAPTIVEGRIDTREE_H_
#define ADAPTIVEGRIDTREE_H_



#include<atomic>
#include<thread>
#include<mutex>
#include<functional>
#include"FastCollisionDetectionLib.h"

// simple dense octree
constexpr int numNodesPerNode = 64;
struct GridTreeNode
{
	GridTreeNode()
	{
		for(int i=0;i<numNodesPerNode;i++)
			childNodes[i]=nullptr;

	}
	std::shared_ptr<GridTreeNode> childNodes[numNodesPerNode];
	std::vector<int> parOrder;

};


template<int concurrency>
class AdaptiveGridTree
{
public:
	AdaptiveGridTree(const float minx, const float miny, const float minz, const float maxx, const float maxy, const float maxz):
		minX(minx),minY(miny),minZ(minz),maxX(maxx),maxY(maxy),maxZ(maxz), workCtr(0)
	{


		for(int i=0;i<concurrency;i++)
		{
			thr.emplace_back([&](){
				FastColDetLib::MemoryPool mem;
				bool work=true;
				while(work)
				{
					auto f = taskQueue.pop();
					if(f)
					{

						f(mem);
						taskCompleteQueue.push2(true);
					}
					else
						break;

				}
				taskCompleteQueue.push(true);
			});


		}
	}

	~AdaptiveGridTree()
	{
		for(int i=0;i<concurrency*2;i++)
			taskQueue.push(nullptr);

		for(int i=0;i<concurrency;i++)
		{
			thr[i].join();
		}
	}

	void clear()
	{
		std::lock_guard<std::mutex> lg(mut);
		parId.reset();
		parMinX.reset();
		parMinY.reset();
		parMinZ.reset();
		parMaxX.reset();
		parMaxY.reset();
		parMaxZ.reset();
	}

	template<typename Derived>
	void addParticles(const int N, Derived * ptr)
	{
		std::lock_guard<std::mutex> lg(mut);
		const int ofsId = parId.allocate(N);
		const int ofsMinx = parMinX.allocate(N);
		const int ofsMiny = parMinY.allocate(N);
		const int ofsMinz = parMinZ.allocate(N);
		const int ofsMaxx = parMaxX.allocate(N);
		const int ofsMaxy = parMaxY.allocate(N);
		const int ofsMaxz = parMaxZ.allocate(N);
		for(int i=0;i<N;i++)
		{
			parId.set(ofsId+i,static_cast<FastColDetLib::IParticle<float>*>(ptr+i)->getId());
			parMinX.set(ofsMinx+i,static_cast<FastColDetLib::IParticle<float>*>(ptr+i)->getMinX());
			parMinY.set(ofsMiny+i,static_cast<FastColDetLib::IParticle<float>*>(ptr+i)->getMinY());
			parMinZ.set(ofsMinz+i,static_cast<FastColDetLib::IParticle<float>*>(ptr+i)->getMinZ());
			parMaxX.set(ofsMaxx+i,static_cast<FastColDetLib::IParticle<float>*>(ptr+i)->getMaxX());
			parMaxY.set(ofsMaxy+i,static_cast<FastColDetLib::IParticle<float>*>(ptr+i)->getMaxY());
			parMaxZ.set(ofsMaxz+i,static_cast<FastColDetLib::IParticle<float>*>(ptr+i)->getMaxZ());
		}
	}


	std::vector<std::pair<int,int>> computeAllPairs(GridTreeNode * parent = nullptr, float minx=0.0f, float miny=0.0f, float minz=0.0f, float maxx=0.0f, float maxy=0.0f, float maxz=0.0f,
			 std::vector<std::pair<int,int>> * result=nullptr, FastColDetLib::Memory<FastColDetLib::FastUnique<int32_t,FastColDetLib::testUniqueLimit>> * mapping = nullptr)
	{

		GridTreeNode root;
		std::vector<std::pair<int,int>> resultRoot;
		if(!parent)
		{

			workCtr=0;
			const int N = parId.size();
			mappingRoot.reset();
			mappingRoot.allocate(N);
			minx = minX;
			miny = minY;
			minz = minZ;
			maxx = maxX;
			maxy = maxY;
			maxz = maxZ;
			parent = &root;



			for(int i=0;i<N;i++)
			{
				parent->parOrder.push_back(i);
				mappingRoot.getRef(i).reset();
			}
			result = &resultRoot;
			mapping = &mappingRoot;

		}



			// build tree until 4096 or less particles left for each leaf
			if(parent->parOrder.size()>10000)
			{
				const int N = parent->parOrder.size();
				for(int i=0;i<numNodesPerNode;i++)
				{

					parent->childNodes[i]=std::make_shared<GridTreeNode>();

				}


				const float stepx = 4.0f/(maxx - minx);
				const float stepy = 4.0f/(maxy - miny);
				const float stepz = 4.0f/(maxz - minz);

				for(int i=0;i<numNodesPerNode;i++)
				{
					accumulator[i].reset();
				}
				for(int i=0;i<N;i++)
				{
					// calculate overlapping region's cell indices
					const int parOrdI = parent->parOrder[i];
					const int mincornerstartx = std::floor((parMinX.get(parOrdI) - minx) * stepx);
					const int maxcornerendx = std::floor((parMaxX.get(parOrdI) - minx) * stepx);
					const int mincornerstarty = std::floor((parMinY.get(parOrdI) - miny) * stepy);
					const int maxcornerendy = std::floor((parMaxY.get(parOrdI) - miny) * stepy);
					const int mincornerstartz = std::floor((parMinZ.get(parOrdI) - minz) * stepz);
					const int maxcornerendz = std::floor((parMaxZ.get(parOrdI) - minz) * stepz);
					for(int ii=mincornerstartz;ii<=maxcornerendz;ii++)
						for(int j=mincornerstarty;j<=maxcornerendy;j++)
							for(int k=mincornerstartx;k<=maxcornerendx;k++)
					{

								if(ii<0 || ii>=4 || j<0 || j>=4 || k<0 || k>=4)
									continue;
								auto & acc = accumulator[k+j*4+ii*16];
								acc.set(acc.allocate(1), parOrdI);
					}
				}

				{

					for(int i=0;i<numNodesPerNode;i++)
					{
						const int n = accumulator[i].size();
						for(int j=0;j<n;j++)
						{
							parent->childNodes[i]->parOrder.push_back(accumulator[i].get(j));
						}
					}
				}

				for(int i=0;i<numNodesPerNode;i++)
				{
					if(parent->childNodes[i]->parOrder.size()>1)
					{
								computeAllPairs(parent->childNodes[i].get(),
										minx+(i&3)/stepx,         miny+((i/4)&3)/stepy,         minz+(i/16)/stepz,
										minx+(i&3)/stepx + 1.0f/stepx, miny+((i/4)&3)/stepy + 1.0f/stepy, minz+(i/16)/stepz + 1.0f/stepz, result,mapping);

					}
				}
			}
			else
			{

				if(parent->parOrder.size()>1)
				{

					// offload to another thread as a sparse-linear-adaptive-grid
					workCtr++;

					taskQueue.push2([&,parent,minx,miny,minz,maxx,maxy,maxz,result](FastColDetLib::MemoryPool  mem)
					{

						FastColDetLib::AdaptiveGridV2 subGrid(mem,minx,miny,minz,maxx,maxy,maxz);
						subGrid.clear();
						{

							subGrid.addParticlesWithoutInterface(parent->parOrder.size(), 0, parent->parOrder,
												parId,
												parMinX, parMinY, parMinZ,
												parMaxX, parMaxY, parMaxZ
												);

							subGrid.buildTree();
						}



						const std::vector<std::pair<int,int>> coll = subGrid.findCollisionsAll();


						if(coll.size()>0)
						{
							std::lock_guard<std::mutex> lg(mut);
							result->insert(result->end(),coll.begin(),coll.end());
						}

					});

				}

			}


			if(parent == &root)
			{
				constexpr int mutN = 1024;
				constexpr int mutN1 = mutN-1;
				FastColDetLib::MutexWithoutFalseSharing mutArr[mutN];

				int endQueue = 0;
				while(endQueue<workCtr)
				{
					if(taskCompleteQueue.pop())
					{
						endQueue++;
					}
				}
				workCtr = 0;
				std::lock_guard<std::mutex> lg(mut);
				resultRootMem.reset();

				{

					{


						const int nr = resultRoot.size();


						for(int i=0;i<nr;i+=1000)
						{
							workCtr++;
							taskQueue.push2([&,i](FastColDetLib::MemoryPool & mem)
							{
								unsigned int lastId = 2000000000; // so you don't have 2 billion particles in a game right?
								mutArr[lastId&mutN1].mut.lock();
								for(int j=0;j<1000;j++)
								{
									if(i+j>=nr)
										break;

									const unsigned int rfirst = resultRoot[i+j].first;
									const int rsecond = resultRoot[i+j].second;
									if(lastId != rfirst)
									{
										//std::lock_guard<std::mutex> lg(mutArr[rfirst&255].mut);
										mutArr[lastId&mutN1].mut.unlock();
										mutArr[rfirst&mutN1].mut.lock();
										mapping->getRef(rfirst).insert(rsecond);
										lastId = rfirst;
									}
									else
									{
										//std::lock_guard<std::mutex> lg(mutArr[rfirst&255].mut);
										mapping->getRef(rfirst).insert(rsecond);
									}
								}

								mutArr[lastId&mutN1].mut.unlock();

							});
						}

						endQueue = 0;
						while(endQueue<workCtr)
						{
							if(taskCompleteQueue.pop())
							{
								endQueue++;
							}
						}
						workCtr = 0;
					}

					resultRoot.clear();
					const int N = mapping->size();

					int allocSize = 0;
					std::vector<std::pair<int,int>> allocOfs;
					for(int i=0;i<N;i++)
					{
						//std::lock_guard<std::mutex> lg(mutArr[i&255].mut);
						const int isz = mapping->getRef(i).size();
						if(isz > 0)
						{
							allocOfs.push_back(std::pair<int,int>(i,allocSize));
							allocSize += isz;
						}
					}

					resultRootMem.allocate(allocSize);

					const int N2 = allocOfs.size();
					for(int i0=0;i0<N2;i0+=10000)
					{
						workCtr++;
						taskQueue.push2([&,i0](FastColDetLib::MemoryPool & mem)
						{

							for(int j=0;j<10000;j++)
							{
								const int i = i0+j;
								if(i>=N2)
									break;


								auto & ref = mapping->getRef(allocOfs[i].first);
								const int n = ref.size();
								const int ofs0 = allocOfs[i].second;
								for(int j=0;j<n;j++)
								{
									resultRootMem.set(ofs0+j,std::pair<int,int>(allocOfs[i].first,ref.get(j)));

								}

							}

						});
					}

					endQueue = 0;
					while(endQueue<workCtr)
					{


						if(taskCompleteQueue.pop())
						{
							endQueue++;
						}

					}
					workCtr = 0;


					resultRoot.resize(resultRootMem.size());
					resultRootMem.writeTo(resultRoot);
				}

			}
			return resultRoot;

	}

private:
	std::mutex mut;
	std::vector<std::thread> thr;
	FastColDetLib::SyncQueue<std::function<void(FastColDetLib::MemoryPool &)>> taskQueue;
	FastColDetLib::SyncQueue<bool> taskCompleteQueue;

	FastColDetLib::Memory<std::pair<int,int>> resultRootMem;
	FastColDetLib::Memory<FastColDetLib::FastUnique<int32_t,FastColDetLib::testUniqueLimit>> mappingRoot;



	const float minX,minY,minZ,maxX,maxY,maxZ;
	int workCtr;
	FastColDetLib::Memory<int> parId;

	FastColDetLib::Memory<float> parMinX;
	FastColDetLib::Memory<float> parMinY;
	FastColDetLib::Memory<float> parMinZ;
	FastColDetLib::Memory<float> parMaxX;
	FastColDetLib::Memory<float> parMaxY;
	FastColDetLib::Memory<float> parMaxZ;


	FastColDetLib::Memory<int> accumulator[numNodesPerNode];

};


#endif /* ADAPTIVEGRIDTREE_H_ */
