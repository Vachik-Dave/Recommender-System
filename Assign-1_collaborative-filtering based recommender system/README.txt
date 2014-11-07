1:- How to:
	compile: make
	Run: ./Assign1 training_filename testing_filename KN_filename output_filename

It takes 4 arguments all of them are filenames.

KN_file format:
	first line: list of K values separated by space.
	second line: list of N values separated by space.
example:
	for K = 3,5,10,20 & N = 5,10,20
/file starts
3 5 10 20
5 10 20
/file ends
same sample is attached in Assignment submission.

Bonus file works similarly with command = ./Assign1Bonus ...

2:-

I have excluded reading/writing file timing & few steps of memory allocation or declaration of varibles.

I calculated similarities only once for all items, so it doesnt changes with K values:
Step 2(a): timing = 222765.00 ms

I used same vector of candidate generation & recommendation list generation for all different N values.
Step 2(b): 
	K = 3	timing = 9973.80 ms
	K = 5	timing = 16945.80 ms
	K = 10 	timing = 31520.10 ms
	k = 20	timing = 53747.80 ms


3:-

I used std sort function from "algorithm" header file in c++.

Here, I've used STL container class, which has vector(linked list) of pair(two tuple).  Each pair is (item_id,similarity value) for every item.

I have written comparison function named "myfunction1 (Line# 280)" which compares simlarity values of pair(two tuple) and return boolean true for first element is greater than second else false.

Using this function and std::sort() we can get sorted list of pairs, where siilarity values are in decreasing order.


4:-

For finding Transpos, I used same approach discussed in lecture. Where First I find count for each column from scanning col_ind vector of CSR. which is O(V) where, V = # nonzero entries in matrix. Then I set row_ptr for Traansposed matrix using above count. Again scanned the val vector of CSR and set corresponding val in Transposed matrix. Complexity would be O(V) where, V = # nonzero entries in matrix OR O(N) where, N = # rows on Transposed matrix. Any of this couls be complexity depending on which is larger.


For finding similarity "get_cos_sim(i,j,Trans)" is used, where cosine similarity is found using same approach discussed in lecture. Here, while finding cosine similarity all nonzero values are 1. Here, I scanned through both item parallelly and check if they are rated by same users, then add them in numerator sum(cosine in code).

The above function is called inside two for loop running from 0 to #items. However, as we know cosine similarity is symmetric I find only one simlarity out of sim(u,v) & sim(v,u) and store both in pair(two tuple) formate described above.
This significantly redurces running time  for large datasets. However theoretically complexity is same i.e. O(N^2) where N = # of items


For storing similarity matrix, Initially I used CSR spaece matrix representation.  Because I stored only half of the Matrix N x N, which is above diagonal as matrix is simmetric for cosine similarity.
However, this method is very slow while retrival(accecting) the elements of the matrix which slower down my preformance in calculating score of Part 2(b). I have attached the screenshot, which shows that part two took almost 50 min for this representation.
This portion is commented in the source code.

Hence, finally I used two dimentional array for N x N matrix for faster access.

5:- 

For different values of K we will get different candidate list, hence from Part 2(b) whole code is in loop of different K values.

For finding candidate list:
for each user:
	- get list of rated items stored in "v"
	- for each rated items:
		- get top K similar items & store it into "candidate" vector
		at the same time check before storing whether the item is already rated or not.
	- remove duplicates from "candidate" vector

For getting Recommendation list:
for each user:
	- get list of rated items stored in "v"
	- for each candidate items:
		- get score: sum through similarity values for rated items (for loop)
			Here, I used already stored similarity matrix "Simi_mat" 2D array. Here accessing elemnets from CSR creats problem.
		- Store score with corresponding candidate item_id as pair(two tuple) in "Reco" vector for each user separately.
	- used same sorting method as above to sort "Reco" vector.


6:- 

I tried inverse of Euclidian Distance, inverse of pearson correlation, and conditional probablity based similarity method (given paper), however most of are same or worse than cosine similarity values.
I've submitted conditional probability based method file.
where before taking Transpose of matrix, all rows are normalized to 1. and all normalized values are take into calculation of similarity.
