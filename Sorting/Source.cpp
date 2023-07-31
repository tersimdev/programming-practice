#include <iostream>
#include <string>

#include "Algorithm.hpp"

using namespace ALGO;


// Author: Terence
// This code demonstrates sorting using insertion, heap, quick, etc.
// Done in C++ using templates.


//MAIN
int main()
{
	//while (1)
	{

		//MIN EDIT DISTANCE
		//std::string inp1, inp2;
		//std::cin >> inp1 >> inp2;
		//std::cout << "Result: " << MinEditDist(inp1, inp2) << std::endl;


		/*SORTING*/

		int arr5[5] = { 3, 2, 4, 1, 5 };
		//int arr20[20] = { 9, 2, 4, 13, 5, 11, 14, 17, 1, 10, 3, 18, 16, 6, 20, 8, 7, 19, 15, 12 };
		float arr20[20] = { 1.2, 1.1, 10, 10, 15, 11, 1.4, 1.7, 1, 10, 3, 1.8, 16, 6, 20, 8, 7, 19, 15, 1.2 };
		//INSERTION SORT
		//InsertionSort(arr, 5);
		//HEAP SORT
		HeapSort(arr20, 20);
		//QUICK SORT
		//QuickSort(arr20, 20);
		//MERGE SORT
		//MergeSort(arr20, 20);
		//HASH SORT
		//HashSort(arr20, 20, 40);
		//for (int i = 0; i < 5; ++i)
			//std::cout << arr5[i] << " ";
		for (int i = 0; i < 20; ++i)
			std::cout << arr20[i] << ", ";
	}
	return 0;
}
