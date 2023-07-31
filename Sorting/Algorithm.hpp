#ifndef ALG_H
#define ALG_H

#include <iostream>
#include <vector>

#define ALGO MyAlgorithms
//Reference: O'REILLY - Algorithms in a Nutshell


namespace MyAlgorithms
{
	//PRIVATE
	//////////////////////////////////////////////////////////////////////////////////////////////
	namespace //anonymous namespace
	{
		//INTERNAL HELPER FUNCTIONS
		//////////////////////////////////////////////////////////////////////////////////////////////
		/***************TABLE OF CONTENTS***************/
		/*CmpFuncMoreThan*/		//Compares if left is more than right, and returns 1, or 0 if same, -1 if less
		/*HashFuncTruncate*/	//Hash function that just truncates the number to nearest int
		/*InternalInsSort*/		//Insertion sort to use internally
		/*Heapify*/				//Helper recursive function for heapsort; validates heap
		/*Partition*/			//Helper partitioning function for quicksort
		/*QuickSortHelper*/		//Helper recursive function for quicksort
		/*MergeSortHelper*/		//Helper recursive function for mergesort
		
		template <typename T>
		int CmpFuncMoreThan(const T left, const T right)
		{
			if (left > right)		return 1;
			else if (right > left)	return -1;
			else					return 0;
		}

		template <typename T>
		int HashFuncTruncate(T value)
		{
			return static_cast<int>(value);
		}

		template <typename T>
		void InternalInsSort(T* arr, const int &arrLength, int(*comparisonFuncPtr)(const T, const T))
		{
			int index;
			T value;
			for (int pos = 1; pos < arrLength; ++pos)
			{
				index = pos - 1;
				value = arr[pos];

				//traverse left side of current value
				while (index >= 0 && comparisonFuncPtr(arr[index], value) == 1)
				{
					arr[index + 1] = arr[index]; //swapping without caring about value of arr[index+1];
					index--; //see next left value
				}
				//insert the value
				arr[index + 1] = value;
			}
		}

		template <typename T>
		void Heapify(T* arr, int(*comparisonFuncPtr)(const T, const T), const int &parentIndex, const int &maxIndex)
		{
			int left = 2 * parentIndex + 1; //index of left child
			int right = 2 * parentIndex + 2; //index of right child
			int largest = 0; //store index of the larger one between parent and two childs

			//find largest
			if (left < maxIndex && comparisonFuncPtr(arr[left], arr[parentIndex]) == 1)
				largest = left;
			else
				largest = parentIndex;
			if (right < maxIndex && comparisonFuncPtr(arr[right], arr[largest]) == 1)
				largest = right;

			if (largest != parentIndex)
			{
				//swap child with parent
				std::swap(arr[parentIndex], arr[largest]);
				//recursively call next heapify with child that was swapped, to validate that branch
				Heapify(arr, comparisonFuncPtr, largest, maxIndex);
			}
		}

		template <typename T>
		int Partition(T* arr, int(*comparisonFuncPtr)(const T, const T), const int &left, const int &right, const int &pivotIndex)
		{
			int store = left;
			T pivotValue = arr[pivotIndex];

			//move pivot to end of array
			std::swap(arr[pivotIndex], arr[right]);

			//move all values <= pivot to left, and all values > pivot to right
			for (int i = left; i < right; ++i)
			{
				if (comparisonFuncPtr(arr[i], pivotValue) == -1)
				{
					std::swap(arr[i], arr[store]);
					store++;
				}
			}
			std::swap(arr[right], arr[store]);

			return store;
		}

		template <typename T>
		void QuickSortHelper(T* arr, int(*comparisonFuncPtr)(const T, const T), const int &left, const int &right, const int &minSize)
		{
			int pivotIndex;
			if (right <= left) 
				return;

			//finding median of three as pivotIndex
			pivotIndex = 1;
			int cmp1, cmp2, cmp3;
			cmp1 = comparisonFuncPtr(arr[0], arr[1]);
			cmp2 = comparisonFuncPtr(arr[1], arr[2]);
			cmp3 = comparisonFuncPtr(arr[0], arr[2]);
			
			if (cmp1 == -1) //2nd more than 1st
			{
				if (cmp2 == 1) //2nd more than 3rd 
				{
					if (cmp3 == 1) //3rd more than 1st 
						pivotIndex = 2;
					else //1st more than 3rd
						pivotIndex = 0;
				}
				else 
					pivotIndex = 1;					
			}
			else //1st more than 2nd
			{
				if (cmp3 == 1) //1st more than 3rd 
				{
					if (cmp2 == 1) //2rd more than 3rd
						pivotIndex = 1;
					else //1st more than 3rd
						pivotIndex = 2;
				}
				else 
					pivotIndex = 0;	
			}

			//partition then re-find index
			pivotIndex = Partition(arr, comparisonFuncPtr, left, right, pivotIndex);

			if (pivotIndex - 1 - left <= minSize)
				InternalInsSort(arr, right + 1, comparisonFuncPtr);
			else
				QuickSortHelper(arr, comparisonFuncPtr, left, pivotIndex - 1, minSize);

			if (right - pivotIndex - 1 <= minSize)
				InternalInsSort(arr, right + 1, comparisonFuncPtr);
			else
				QuickSortHelper(arr, comparisonFuncPtr, pivotIndex + 1, right, minSize);
		}
		
		template <typename T>
		void MergeSortHelper(T* arr, T* result, int start, int end, int(*comparisonFuncPtr)(const T, const T))
		{
			if (end - start < 2)
				return; //escape case
			else if (end - start == 2)
			{
				//sort result
				if (comparisonFuncPtr(result[start], result[start + 1]) == 1)				
					std::swap(result[start], result[start + 1]);
				return;					
			}

			int mid = (start + end) / 2;

			MergeSortHelper(result, arr, start, mid, comparisonFuncPtr);
			MergeSortHelper(result, arr, mid, end, comparisonFuncPtr);

			//merge left and right side
			int i = start;
			int j = mid;
			for (int idx = start; idx < end; ++idx)
			{
				if (j >= end || (i < mid && arr[i] < arr[j]))
				{
					result[idx] = arr[i];
					i++;
				}
				else
				{
					result[idx] = arr[j];
					j++;
				}
			}

		}

		//UTILITY STATELESS FUNCTIONS
		/*Everything below here should have NO state*/
		//////////////////////////////////////////////////////////////////////////////////////////////
		/***************TABLE OF CONTENTS***************/
		/*MathExt*/			//Math functions
		/*Mtx2D*/			//Helper functions for a 2D Matrix using double pointer

		namespace MathExt
		{
			template <typename T>
			inline static T Min(const T &t1, const T &t2)
			{
				if (t1 <= t2)
					return t1;
				else
					return t2;
			}

		}

		class Mtx2D
		{
		public:
			template <typename T>
			inline static void Create(T** &arr, int w, int h)
			{
				arr = new T*[w];

				for (int i = 0; i < w; ++i)
				{
					arr[i] = new T[h];
					for (int j = 0; j < h; ++j)
					{
						arr[i][j] = 0;
					}
				}
			}

			template <typename T>
			inline static void Print(T** arr, int w, int h)
			{
				for (int i = 0; i < w; ++i)
				{
					for (int j = 0; j < h; ++j)
					{
						std::cout << arr[i][j] << " ";
					}
					std::cout << std::endl;
				}
				std::cout << std::endl;
			}
		};
	}

	
	//PUBLIC API
	//////////////////////////////////////////////////////////////////////////////////////////////
	/***************TABLE OF CONTENTS***************/
	/*MinEditDist*/			//Calculates minimum amount of edits to make two strings the same
	/*InsertionSort*/		//Use to sort a small array (<10) or almost sorted array
	/*HeapSort*/			//Use to avoid worst-case of QuickSort; Generally use this or QuickSort
	/*QuickSort*/			//Use for good average-case sort; Generally use this or HeapSort
	/*HashSort*/			//Use when elements are dense and uniformly distributed
	/*MergeSort*/			//Use when stable sort is required

	int MinEditDist(const std::string &s1, const std::string &s2)
	{
		//lengths of string including the empty "" string
		int w = s1.length() + 1;
		int h = s2.length() + 1;

		//create two-d matrix m
		int** m;
		Mtx2D::Create(m, w, h);

		//init m
		for (int i = 0; i < w; ++i)
			m[i][0] = i;
		for (int j = 0; j < h; ++j)
			m[0][j] = j;

		//calculate cost for each node
		for (int i = 1; i < w; ++i)
		{
			for (int j = 1; j < h; ++j)
			{
				//calc costs
				int replaceCost = m[i - 1][j - 1] + ((s1[i - 1] == s2[j - 1]) ? 0 : 1);
				int removeCost = m[i - 1][j] + 1;
				int insertCost = m[i][j - 1] + 1;

				//finding minimum
				m[i][j] = MathExt::Min(
					MathExt::Min(replaceCost, removeCost),
					MathExt::Min(replaceCost, insertCost)
				);
			}
		}

		Mtx2D::Print(m, w, h); //for debugging

		int ret = m[w - 1][h - 1];
		delete[] m;
		return ret;
	}

	template<typename T>
	void InsertionSort(T* arr, const int &arrLength, int(*comparisonFuncPtr)(const T, const T) = &CmpFuncMoreThan)
	{
		InternalInsSort(arr, arrLength, comparisonFuncPtr);
	}

	template<typename T>
	void HeapSort(T* arr, const int &arrLength, int(*comparisonFuncPtr)(const T, const T) = &CmpFuncMoreThan)
	{
		//build heap, a binary tree where parent > child
		for (int i = arrLength * 0.5 - 1; i >= 0; --i)
		{
			//starting from lowest parent, build and validify branch in-place
			Heapify(arr, comparisonFuncPtr, i, arrLength);
		}

		//sort root(biggest number) to the right of the array
		for (int i = arrLength - 1; i >= 1; --i)
		{
			//swap root with bottom-most leaf
			std::swap(arr[0], arr[i]);
			//re-evaluate heap, ignoring bottom-most leaf
			Heapify(arr, comparisonFuncPtr, 0, i);
		}
	}

	template<typename T>
	void QuickSort(T* arr, const int &arrLength, int(*comparisonFuncPtr)(const T, const T) = &CmpFuncMoreThan)
	{
		QuickSortHelper(arr, comparisonFuncPtr, 0, arrLength - 1, 7);
	}

	template<typename T>
	void HashSort(T* arr, const int &arrLength, const int &numBuckets, int(*comparisonFuncPtr)(const T, const T) = &CmpFuncMoreThan, int(*hashFuncPtr)(const T) = &HashFuncTruncate)
	{
		//create buckets
		T** buckets;
		int maxBucketSize = arrLength / numBuckets * 2;
		Mtx2D::Create(buckets, numBuckets, maxBucketSize);
		int bucketId = 0; //current bucket index
		int arrayId = 0; //curent array index for extraction
		int* bucketSizes = new int[numBuckets]; //current bucket (filled) size

		//init bucket sizes
		for (int i = 0; i < numBuckets; ++i)
		{
			bucketSizes[i] = 0;
		}
		
		//assign data to buckets
		for (int i = 0; i < arrLength; ++i)
		{
			bucketId = hashFuncPtr(arr[i]);
			buckets[bucketId][bucketSizes[bucketId]++] = arr[i]; //basically pushing back
			if (bucketSizes[bucketId] >= 10)
				std::cout << "Too much data in bucket " << bucketId << std::endl;
		}

		//sort buckets and extract
		for (bucketId = 0; bucketId < numBuckets; ++bucketId)
		{
			InternalInsSort(buckets[bucketId], bucketSizes[bucketId], comparisonFuncPtr);

			//extract
			for (int i = 0; i < bucketSizes[bucketId]; ++i)
			{
				arr[arrayId++] = buckets[bucketId][i];
			}
		}

		delete bucketSizes;
		delete[] buckets;
	}

	template<typename T>
	void MergeSort(T* arr, const int &arrLength, int(*comparisonFuncPtr)(const T, const T) = &CmpFuncMoreThan)
	{
		//create copy of array
		T* copy = new T[arrLength];
		for (int i = 0; i < arrLength; ++i)
			copy[i] = arr[i];

		MergeSortHelper(copy, arr, 0, arrLength, comparisonFuncPtr);
	}


}
#endif