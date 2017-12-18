//
// (c) MERL 2012 - 2013
//
/* 	
	Utility functions for sorting and finding the minimum
*/

#pragma once

#include <vector>
#include <algorithm>

template <typename T>
class IndexSorter
{
private:
	struct Pair
	{
		T val;
		int idx;
	};

	static bool Compare(const Pair& lhs, const Pair& rhs)
	{
		return (lhs.val < rhs.val);
	};

	static bool decCompare(const Pair& lhs, const Pair& rhs)
	{
		return (lhs.val > rhs.val);
	};

public:
	static void Sort(T* valVec, const int n, int* idxVec)
	{
		std::vector<Pair> pairVec(n);
		for (int i=0; i<n; ++i)
		{
			pairVec[i].val = valVec[i];
			pairVec[i].idx = idxVec[i];
		}

		std::sort(pairVec.begin(), pairVec.end(), Compare);

		for (int i=0; i<n; ++i)
		{
			valVec[i] = pairVec[i].val;
			idxVec[i] = pairVec[i].idx;
		}
	};

	static void decSort(T* valVec, const int n, int* idxVec)
	{
		std::vector<Pair> pairVec(n);
		for (int i=0; i<n; ++i)
		{
			pairVec[i].val = valVec[i];
			pairVec[i].idx = idxVec[i];
		}

		std::sort(pairVec.begin(), pairVec.end(), decCompare);

		for (int i=0; i<n; ++i)
		{
			valVec[i] = pairVec[i].val;
			idxVec[i] = pairVec[i].idx;
		}
	};

	static const int MinIndex(T* valVec, const int n)
	{
		int idx = 0;
		T minVal = valVec[0];

		for (int i=1; i<n; ++i)
		{
			if (valVec[i] < minVal)
			{
				minVal = valVec[i];
				idx = i;
			}
		}

		return idx;
	};
};


