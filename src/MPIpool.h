#pragma once
#include <thread>
#include "MPIwrapper.h"

template<typename TaskType, typename ResType>
class MPIpool {
public:
	MPIpool(MPI &mpi, const vector<TaskType>& task, vector<ResType> &result, const size_t chunkSize) : mpi(mpi) {
		if (mpi.gridSize < 2) throw runtime_error("You must run at least 2 MPI tasks.");
		if (!mpi.isRoot()) return;
		result.resize(task.size());
		poolMaster(mpi, task, result, chunkSize);
	}

	vector<TaskType> getTask() {
		vector<TaskType> tmp;
		mpi.recv(tmp, mpi.root);
		return tmp;
	}
	void submit(const vector<ResType>& res) {
		mpi.send(res, mpi.root);
	}

private:
	MPI &mpi;

	static void poolMaster(MPI &mpi, const vector<TaskType>& task, vector<ResType>& result, const size_t chunkSize) {
		size_t curIt = 0;
		std::mutex m;
		Stopwatch tmr;
		auto nextRange = [&]() -> pair<size_t, size_t> {
			m.lock();
			const size_t from = curIt;
			const size_t to = (curIt += min(chunkSize, task.size() - curIt));
			m.unlock();
			
			const size_t tn = from / chunkSize;
			const size_t tt = task.size() / chunkSize;
			if(tn < tt) cout << "Task dispatched " << tn+1 << " / " << tt << " : " << tmr.stop() << "sec." << endl;
			else cout << "Worker relieved." << endl;

			return { from, to };
		};

		auto ker = [&](const int targetId) {
			while (1) {
				auto r = nextRange();
				const size_t chunkSz = r.second - r.first;
				mpi.send(task.data() + r.first, chunkSz, targetId);
				if (!chunkSz) break;
				mpi.recv(result.data() + r.first, chunkSz, targetId);
			}
		};

		vector<std::thread> workers(mpi.gridSize - 1);
		for (int i = 0; i < workers.size(); ++i)
			workers[i] = std::thread(ker, i + 1);
		for (int i = 0; i < workers.size(); ++i)
			workers[i].join();

	}
};