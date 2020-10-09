#pragma once

#include <mpi.h>
#include <exception>
#include <sstream>
#include <vector>
#include <string>
#include <cstring>
#include <climits>

#define HOSTNAME_BUFF_SZ 50
class MPIwrapper {
public:
	int myId;
	int gridSize;
	const int root = 0;
	std::string host;

	struct nodeInfo {
		int id;
		char hostname[HOSTNAME_BUFF_SZ];
	};
	nodeInfo me;

	MPIwrapper(MPI_Comm comm = MPI_COMM_WORLD) : comm(comm) {
		int pr;
		MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &pr);
		if (pr != MPI_THREAD_MULTIPLE) throw std::runtime_error("Multithread MPI request failed.");
		MPI_Comm_size(comm, &gridSize);
		MPI_Comm_rank(comm, &myId);

		//get host name
		char hostname[200];
		int len;
		MPI_Get_processor_name(hostname, &len);
		host = std::string(hostname);

		me.id = myId;
		strncpy(me.hostname, hostname, HOSTNAME_BUFF_SZ-1);
		me.hostname[HOSTNAME_BUFF_SZ-1] = '\0';
	}

	virtual ~MPIwrapper() {
		MPI_Finalize();
	}

	class Error : public std::exception {
	public:
		Error(const int rank, const int e) {
			char error_string[BUFSIZ];
			int length_of_error_string, error_class;
			MPI_Error_class(e, &error_class);
			MPI_Error_string(error_class, error_string, &length_of_error_string);
			std::stringstream ss;
			ss << "MPI error: rank " << rank << ", err class: " << error_string;
			MPI_Error_string(e, error_string, &length_of_error_string);
			ss << ", err: " << error_string;
			error = ss.str();
		}
		virtual const char* what() const noexcept {
			return error.c_str();
		}
	private:
		std::string error;
	};

	std::tuple<int, bool> localId() {
		static int localId = 0;
		static bool rootMachine = false;
		static bool done = false;
		if(done) return std::make_tuple(localId, rootMachine);
		done = true;

		char rootHostname[HOSTNAME_BUFF_SZ];
		strcpy(rootHostname, me.hostname);
		Bcast(rootHostname);
		rootMachine = strcmp(rootHostname, me.hostname) == 0;

		std::vector<nodeInfo> nodes;
		Allgather(me, nodes);
		std::vector<nodeInfo> locals;
		std::copy_if(nodes.begin(), nodes.end(), back_inserter(locals), [&](const nodeInfo &n) {return strcmp(n.hostname, me.hostname) == 0; });
		localId = std::find_if(locals.begin(), locals.end(), [&](const nodeInfo &n) {return n.id == me.id; }) - locals.begin();
		return std::make_tuple(localId, rootMachine);
	}

	char *hostname() {
		static char hostname[100];
		int len;
		MPI_Get_processor_name(hostname, &len);
		return hostname;
	}

	template <typename T>
	void Bcast(std::vector<T> &v) {
		size_t sz = v.size();
		Bcast(sz);
		v.resize(sz);
		error(MPI_Bcast64((uint8_t*)v.data(), sz * sizeof(T), MPI_UINT8_T, root, comm));
	}

	template <typename T>
	void Bcast(T &v) {
		error(MPI_Bcast(&v, sizeof(T), MPI_UINT8_T, root, comm));
	}
	template <typename T>
	void Bcast(T &&v) {
		Bcast(v);
	}

	template<typename T, typename... Args>
	void Bcast(T &v, Args&... args) {
		Bcast(v);
		Bcast(args...);
	}

	template <typename T>
	void Allgather(const T &v, std::vector<T> &res) {
		res.resize(gridSize);
		error(MPI_Allgather(&v, sizeof(T), MPI_UINT8_T, res.data(), sizeof(T), MPI_UINT8_T, comm));
	}

	template <typename T>
	void Scatter(std::vector<T> &src, std::vector<T> &part, const int sz[]) {
		part.resize(sz[myId]);
		std::vector<int> sizes(sz, sz + gridSize);
		for (auto &s : sizes) s *= sizeof(T);
		std::vector<int> displacement(sizes.size());
		for (size_t i = 0; i < sizes.size() - 1; ++i)
			displacement[i + 1] = displacement[i] + sizes[i];
		error(MPI_Scatterv(src.data(), sizes.data(), displacement.data(), MPI_UINT8_T, part.data(), sizes[myId], MPI_UINT8_T, root, comm));
	}

	template <typename T>
	void Gather(std::vector<T> &v, std::vector<T> &res) {
		std::vector<int> sizes;
		Allgather(int(v.size()), sizes);
		for (auto &s : sizes) s *= sizeof(T);
		std::vector<int> displacement(sizes.size());
		for (size_t i = 0; i < sizes.size() - 1; ++i)
			displacement[i + 1] = displacement[i] + sizes[i];
		const size_t totalSize = (*displacement.crbegin() + *sizes.crbegin()) / sizeof(T);
		res.resize(totalSize);
		error(MPI_Gatherv(v.data(), sizes[myId], MPI_UINT8_T, res.data(), sizes.data(), displacement.data(), MPI_UINT8_T, root, comm));
	}

	void Barrier() {
		MPI_Barrier(comm);
	}
	
	template <typename T>
	void send(const T* v, const size_t sz, const int destId) {
		send(sz, destId);
		if (!sz) return;
		error(MPI_Send(v, sz*sizeof(T), MPI_UINT8_T, destId, 0, comm));
	}
	
	template <typename T>
	void recv(T* v, const size_t sz, const int srcId) {
		size_t rsz;
		recv(rsz, srcId);
		if (rsz != sz) throw std::runtime_error("MPI::recv(T* v, const size_t sz, const int srcId) - size of v differ");
		if (!sz) return;
		error(MPI_Recv(v, sz*sizeof(T), MPI_UINT8_T, srcId, 0, comm, MPI_STATUS_IGNORE));
	}
	
	template <typename T>
	void send(const std::vector<T>& v, const int destId) {
		const size_t sz = v.size();
		send(sz, destId);
		if (!sz) return;
		error(MPI_Send(v.data(), sz*sizeof(T), MPI_UINT8_T, destId, 0, comm));
	}
	
	template <typename T>
	void recv(std::vector<T>& v, const int srcId) {
		size_t sz;
		recv(sz, srcId);
		v.resize(sz);
		if (!sz) return;
		error(MPI_Recv(v.data(), sz*sizeof(T), MPI_UINT8_T, srcId, 0, comm, MPI_STATUS_IGNORE));
	}	

	template <typename T>
	void send(const T& v, const int destId) {
		error(MPI_Send(&v, sizeof(T), MPI_UINT8_T, destId, 0, comm));
	}
	
	template <typename T>
	void recv(T& v, const int srcId) {
		error(MPI_Recv(&v, sizeof(T), MPI_UINT8_T, srcId, 0, comm, MPI_STATUS_IGNORE));
	}
	
	bool isRoot() { return myId == root; }
	bool isLocalRoot() { return std::get<0>(localId()) == 0; }

private:
#define MPI_MAX_BCAST size_t(INT_MAX/4)
	template <typename T>
	int MPI_Bcast64(T* buffer, size_t count, MPI_Datatype datatype, int root, MPI_Comm comm) {
		for (; count > MPI_MAX_BCAST; count -= MPI_MAX_BCAST, buffer += MPI_MAX_BCAST)
			if (const int e = MPI_Bcast(buffer, MPI_MAX_BCAST, MPI_UINT8_T, root, comm) != MPI_SUCCESS) return e;
		return MPI_Bcast(buffer, count, MPI_UINT8_T, root, comm);
	}

	MPI_Comm comm;
	void error(const int e) {
		if (e != MPI_SUCCESS)
			throw Error(myId, e);
	}
};
