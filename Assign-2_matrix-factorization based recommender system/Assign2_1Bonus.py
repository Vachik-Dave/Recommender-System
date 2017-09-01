import os;
import sys;
import numpy as np;
import datetime as dt;
import math;
#from numpy import linalg as LA;


class CSR:
	# Different ways of initializations
	def	__init__(self):
		self.R = 0;
		self.C = 0;
		self.V = 0;

	def	File_init(self,filetoread):
		self.file_id = filetoread;

	def	My_init(self,a,b,c):
		self.R = a;
		self.C = b;
		self.V = c;
		self.val = [0] * self.V;
		self.col_ind = [0] * self.V;
		self.row_ptr = [0] * (self.R+1);
	# Get specific element
	def	Get(r,c):
		for i in range(self.row_ptr[r],self.row_ptr[r+1]):
			if self.col_ind[i] == c:
				return self.val[i];
			if self.col_ind[i] > c:
				return 0.0;
		return 0.0;		

	# Get specific element
	def	Get_mean(self):
		sum1 = 0.0;
		Max = 0.0;
		Min = 99999.0;
		for i in range(0,self.V):
			sum1 = sum1 + self.val[i]; 
			if Max < self.val[i]:
				Max = self.val[i];
			if Min > self.val[i]:
				Min = self.val[i];
		self.mean = sum1 / self.V;
		return [self.mean,Max,Min];
	def	Get_allMean(self):
		Max = [0.0]*self.R;
		Min = [99999.0]*self.R;
		self.mean = 0.0;
		self.userR = [];
		self.itemR = [0.0] * self.C;
		col_count = [0.0] * self.C;
		for 	i 	in 	range(0,self.R):
			sum1 = 0.0;
			count = self.row_ptr[i+1]-self.row_ptr[i];
			for j in range(self.row_ptr[i],self.row_ptr[i+1]):
				sum1 = sum1 + self.val[j]; 
				c = self.col_ind[j];
				self.itemR[c] = self.itemR[c] + self.val[j];
				col_count[c] = col_count[c]+1;
				self.mean = self.mean + self.val[j];
				if Max[i] < self.val[j]:
					Max[i] = self.val[j];
				if Min[i] > self.val[j]:
					Min[i] = self.val[j];
			if count != 0:
				self.userR.append(sum1 / count);
		
		for i in range(0,len(self.itemR)):
			if col_count[i] != 0:
				self.itemR[i] = self.itemR[i] / col_count[i]; 
		self.mean = self.mean / self.V;
		return [self.mean,self.userR,self.itemR,Max,Min];

	def	Adjust_by_mean(self):
		for i in range(0,self.V):
			self.val[i] = self.val[i] - self.mean;
	def	Adjust_Bias(self,mul):
		for 	i 	in 	range(0,self.R):
			for j in range(self.row_ptr[i],self.row_ptr[i+1]):
				col = self.col_ind[j];
				bias = self.userR[i] + self.itemR[col];
				if self.itemR[col] == 0:
					bias = self.userR[i] + self.mean;
				self.val[j] = self.val[j] - (mul*bias); #- self.mean;
		
	# Read file and assign data and class values
	def	readFile(self):
		# stores each instance of the database
		self.val = [];
		self.col_ind = [];
		self.row_ptr = [];
		self.row_ptr.append(0);
		count =0;
		for	line	in	self.file_id:
			tuplex= line.strip().split(" ");
			if self.R == 0:
				self.R = int(tuplex[0]);
				self.C = int(tuplex[1]);
				self.V = int(tuplex[2]);
				continue;
			for	position	in	range(0,len(tuplex)):
				if position % 2 == 0:
					self.col_ind.append(int(tuplex[position])-1);
				else:
					self.val.append(int(tuplex[position]));
					count = count +1;
			self.row_ptr.append(count);

	#display CSR data
	def	display(self):
		print "#Rows = ";
		print self.R;
		print "\n#Cols = ";
		print self.C;
		print "\n#NonZero = ";
		print self.V;
		print "\nCol_ind: ";
		print self.col_ind;
		print "\nValues: ";
		print self.val;
		print "\nRow_ptr:";
		print self.row_ptr;

#End of CSR class
#############################################################################################################

def Transpose(M):

	Trans = CSR();	
	Trans.My_init(M.C,M.R,M.V);
	row_counts2 = [0] * Trans.R;
	
	for i  in range(0,Trans.V):
		t = M.col_ind[i];
		row_counts2[t]+= 1;

	for i in range(1,Trans.R+1):
		Trans.row_ptr[i] = Trans.row_ptr[i-1] + row_counts2[i-1];

	row_counts2 = [0] * Trans.R;
	for i in range(0,M.R):
		for j in range(M.row_ptr[i],M.row_ptr[i+1]):
			c = M.col_ind[j]; 
			if row_counts2[c] > (Trans.row_ptr[c+1]-Trans.row_ptr[c]):
				print "Error: count problem in trnspose\n";
			Trans.col_ind[Trans.row_ptr[c] + row_counts2[c]] = i; 
			Trans.val[Trans.row_ptr[c] + row_counts2[c]] = M.val[j];  
			row_counts2[c] += 1;  
	return Trans;
#############################################################################################################

def LS_closed(X,Y,K,Lemda):
	result = [];
	outer_list = [];
	for i in range(0,Y.shape[0]):
		t = Y[i][:];
		mul = np.outer(t,t);
		outer_list.append(mul);
	for i in range(0,X.R):
		st = X.row_ptr[i];
		end = X.row_ptr[i+1];
		sum1 = np.zeros((1,K), dtype=np.float);
		sum2 = np.zeros((K,K), dtype=np.float);
		ret = np.array([]);
		for j in range(st,end):
			b = X.col_ind[j];
#			value = X.val[j];
#			print b;
#			print Y.shape

			s1 = Y[b][:] * X.val[j];
			sum1 = sum1 + s1;
			
#			s2 = np.outer(t,t);
			s2 = outer_list[b];
			sum2 = sum2 + s2;
		I = np.matrix(np.identity(sum2.shape[0]));
		I =  I * Lemda;
		r = I + sum2;
#		print r;
		secondT = np.linalg.inv(r);
#		print secondT;
		res = np.dot(sum1,secondT);
#		print res;
		result.append(res);

	ret = np.vstack(result);
#	print ret.shape
	return ret;

#############################################################################################################

def Get_objVal(Mat,P,Q,Lemda):
	norm1 = np.linalg.norm(P, ord=2);
	normP = norm1**2;
	norm2 = np.linalg.norm(Q, ord=2);
	normQ = norm2**2;
	sum1 = 0.0;	
	for i in range(0,Mat.R):
		st = Mat.row_ptr[i];
		end = Mat.row_ptr[i+1];
		pi = P[i][:];
		for j in range(st,end):
			b = Mat.col_ind[j];
			qj = Q[b][:];
			v = np.dot(pi,qj.T);
			res = (Mat.val[j]-v)**2;
			sum1 = sum1 + res;
	ret = sum1 + Lemda*(normP+normQ);
	ret = ret /2.0;
	return ret;



#############################################################################################################

def Get_MSE(Test,P,Q,mean,userR,itemR,Bmul,Max,Min):
	MSE = 0.0;
	for i in range(0,Test.R):
		st = Test.row_ptr[i];
		end = Test.row_ptr[i+1];
		pi = P[i][:];
		for j in range(st,end):
			b = Test.col_ind[j];
			qj = Q[b][:];
			Pred = np.dot(pi,qj.T);
			bias =  userR[i] + itemR[b];
			if  itemR[b] == 0:
				bias = userR[i]+mean;
			Pred = Pred + (Bmul*bias);#+ mean;
			if Pred > Max[i]:
				Pred = Max[i];
			if Pred < Min[i]:
				Pred = Min[i];
			SE = (Test.val[j]-Pred)**2;
			MSE = MSE + SE;
	MSE = MSE / Test.V;

	return MSE;
#############################################################################################################
def main():

#	parser = argparse.ArgumentParser();
#	parser.add_argument("-i","--inputfile");
#	parser.add_argument("-t","--testingfile");
#	parser.add_argument("-p","--parameterfile");
#	parser.add_argument("-o","--outputfile");
#	args = parser.parse_args();	
 
	if	len(sys.argv) < 4:
		print 'Please give the filenames for training, testing, parameters'+os.linesep;
		sys.exit(1);
	try:	
		#file name
		train_name = sys.argv[1];
		train_file = open(train_name);
		test_name = sys.argv[2];
		test_file = open(test_name);
		para_name = sys.argv[3];
		para_file = open(para_name);
		out_file = open('result.txt');
	except IOError,IndexError:
		print	'File cannot open OR Bad file name'+os.linesep;
		sys.exit(1);

	#minsup = 2;

	# CSR object
	Mat = CSR();
	Mat.File_init(train_file);
	#reading the file
	Mat.readFile();
	#display
#	Mat.display();

	[mean,userR,itemR,Max,Min] = Mat.Get_allMean();
	#bias_multiplier 

#	Bmul = 0.05;
	Bmul = 0.5;

	Mat.Adjust_Bias(Bmul);
#	Mat.Adjust_by_mean();

	# Get Transposed matrix
	Trans = Transpose(Mat);
#	Trans.display();


	# CSR object
	Test = CSR();
	Test.File_init(test_file);
	#reading the file
	Test.readFile();

	K_list = [];
	L_list = [];
	Iter_list = [];
	eps_list = [];
	t = 0;
	for	line	in	para_file:
		tuplex= line.strip().split(" ");
		for	pos	in	range(0,len(tuplex)):
			if t == 0:
				K_list.append(int(tuplex[pos]));
			if t == 1:
				L_list.append(float(tuplex[pos]))
			if t == 2:
				Iter_list.append(int(tuplex[pos]));
			if t == 3:
				eps_list.append(float(tuplex[pos]))
		t = t +1;


#	K = 20;
#	Lemda = 3;
#	maxIter = 100;
#	epsilon = 0.00001;



	if os.path.isfile('Bonus_result.txt'):
		out_file = open('Bonus_result.txt');
		result = [];
		for 	line 	in 	out_file:
			result.append(line);
		out_file = open('Bonus_result.txt','w');
		for 	l 	in 	result:
			out_file.write(l);

	else:
		out_file = open('Bonus_result.txt','w');

#	out_file.write('\n');

	for e in range(0,len(K_list)):	
		K = K_list[e];
	
		for f in range(0,len(L_list)):
			Lemda = L_list[f];
			for g in range(0,len(Iter_list)):
				maxIter = Iter_list[g];
				for h in range(0,len(eps_list)):
					epsilon = eps_list[h];

# Repeat the whole process for all combinations
##########################################################################

					P = np.empty((Mat.R,K));
					P[:] = 1.0/K;
					Q = np.empty((Mat.C,K));
					Q[:] = 1.0/K;

					t = 0;
					prev = Get_objVal(Mat,P,Q,Lemda);

					# Start the timmer
					start_time = dt.datetime.now();		

					min_MSE = 99999.0;
					while t < maxIter:
						t = t +1;
						P = LS_closed(Mat,Q,K,Lemda);
						Q = LS_closed(Trans,P,K,Lemda);
						cur = Get_objVal(Mat,P,Q,Lemda);
		
						check = abs(cur-prev);
						check = check/prev;
						if check < epsilon:
							break;
						M = Get_MSE(Test,P,Q,mean,userR,itemR,Bmul,Max,Min);
						if M < min_MSE:
							min_MSE = M;
						print  str(t) + " -> "+ str(cur) + " - " + str(M);  
					end_time = dt.datetime.now();
					time = end_time - start_time;
					print 'Training time = '+ str(time.total_seconds()) + ' Seconds';

					# Start the timmer
#					start_time = dt.datetime.now();		

#					MSE = Get_MSE(Test,P,Q,meanR,Max,Min);
#					print MSE;
					RMSE = math.sqrt(min_MSE);
#					print RMSE;

#					end_time = dt.datetime.now();
#					time = end_time - start_time;
#					print 'Testing time = '+ str(time.total_seconds()) + ' Seconds';
					tm = time.total_seconds() * 1000.00;
					time_msec = float("{0:.2f}".format(tm));
					print repr(K) + ' ' + repr(Lemda)+ ' '+ repr(maxIter) +' '+repr(epsilon)+ ' ' +str(float(min_MSE))+ ' '+ repr(RMSE)+' '+str(time_msec)+' secs\n';
					out_file.write(repr(K) + ' ' + repr(Lemda)+ ' '+ repr(maxIter) +' '+repr(epsilon)+ ' ' +str(float(min_MSE))+ ' '+ repr(RMSE)+' '+str(time_msec)+'\n');
#End of main function
###################################################################################################################
if	__name__== "__main__":
	main();
