#pragma once

template <typename PointType, typename DataType>
class DataFile {
public:
	struct Sample {
		PointType p;
		DataType v;
	};
private:
	using A = std::allocator<T>;

public:
	typedef A allocator_type;
	typedef typename A::value_type value_type;
	typedef typename A::reference reference;
	typedef typename A::const_reference const_reference;
	typedef typename A::difference_type difference_type;
	typedef typename A::size_type size_type;

	class iterator {
	public:
		typedef typename A::difference_type difference_type;
		typedef typename A::value_type value_type;
		typedef typename A::reference reference;
		typedef typename A::pointer pointer;
		typedef std::random_access_iterator_tag iterator_category;

		iterator() {}
		//iterator(const iterator&);
		virtual ~iterator() {}
		//iterator& operator=(const iterator&);
		virtual bool operator==(const iterator&) const {}
		bool operator!=(const iterator& it) const { return !(*this == it); }
		virtual iterator& operator++() {}
		virtual reference operator*() const { return reference(); }
		virtual pointer operator->() const { return pointer(); }

	};
};