#pragma once

#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
using namespace boost::interprocess;

template <class T>
using ShmemAllocator = boost::interprocess::allocator<T, managed_shared_memory::segment_manager>;

template <class VT>
class SharedMemBase {
public:
	virtual ~SharedMemBase() {}
	VT* data;
};

template <class VT>
class SharedMemMaster : public SharedMemBase<VT> {
	struct shm_remove
	{
		shm_remove() { shared_memory_object::remove("hexElSharedMemory"); }
		~shm_remove() { shared_memory_object::remove("hexElSharedMemory"); }
	} remover;
	managed_shared_memory segment;
	const ShmemAllocator<typename VT::value_type> alloc_inst;

public:
	SharedMemMaster(const size_t sz) :
		segment(create_only, "hexElSharedMemory", sz * sizeof(typename VT::value_type) + 0x10000),
		alloc_inst(segment.get_segment_manager()) {
		SharedMemBase<VT>::data = segment.construct<VT>("hexElVector")(alloc_inst);
	}

	~SharedMemMaster() {
		segment.destroy<VT>("hexElVector");
	}
};

template <class VT>
class SharedMemSlave : public SharedMemBase<VT> {
	managed_shared_memory segment;
public:
	SharedMemSlave() :
		segment(open_only, "hexElSharedMemory") {
		SharedMemBase<VT>::data = segment.find<VT>("hexElVector").first;
		cout << "SLAVE SHARED MEMORY CREATED!" << endl;
	}
};