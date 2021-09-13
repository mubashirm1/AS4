from math import sqrt
#Reading a matrix from a text file and returns it as a list
def read_Matrix(file_name):
    file=open(file_name,'r')
    matrix=[]
    for i in file:
        row_txt=i.split()
        row = [float(j) for j in row_txt]
        matrix.append(row)
    file.close
    return(matrix)

#Writes the matrix to a text file
def write_Matrix(file_name,output):
    file=open(file_name,'w')
    mat=[]
    for i in output:
        row = ['{:.2f}'.format(j) for j in i]
        mat.append('\t'.join(row))
    file.writelines('\n'.join(mat))
    file.close

#Displays the matrix
def display_Matrix(output):
    for i in output:
        row = ['{:.2f}'.format(j) for j in i]
        print('\t'.join(row))

#Coverting matrix into an augmented matrix
def aug_Matrix(a,b):
    for i in range(len(a)):
        for j in range(len(b[i])):
            a[i].append(b[i][j])
    return a


#Extracting the result from augmented matrix
def extract_Matrix(ab,b):
    a=[]
    column_no=len(ab[0])-len(b[0])
    for i in range(column_no):
        b[i]=ab[i][column_no:]
        a[i]=ab[i][:column_no]
    return a,b


#Partial pivot
def partial_Pivot(ab,m):
    counter=0
    if(ab[m][m]==0 or abs(ab[m][m]) < 10 ** (-12)):#Checks if the value at pivot point is zero or close to zero
        max=m
        f=0
        for i in range(m+1,len(ab)):#Finds the greatest value below the pivot point in the column 
            if(ab[i][m]>max):
                max=i
                f=1
        if(f!=1):
            max=len(ab)-1
        temp=ab[max]#swapping the rows
        ab[max]=ab[m]
        ab[m]=temp
        counter=1
    return ab,counter 


#Guass Jordan meathod
def gauss_Jordan(a,b):
    ab=aug_Matrix(a,b)
    n_columns=len(ab[0])
    n_rows=len(ab)
    for i in range(n_rows):
        ab,_=partial_Pivot(ab,i)
        pivot=ab[i][i]
        for j in range(n_columns):
            ab[i][j]=ab[i][j]/pivot
        pivot=ab[i][i]
        for k in range(n_rows):
            if k!=i:
                sub=ab[k][i]
                for l in range(n_columns):
                    ab[k][l]-=(sub/pivot)*ab[i][l]
    _,result=extract_Matrix(ab,b)
    return(result)    


#Method to find the determinant
def determinant(a):
    n_rows=len(a)
    count=0
    for i in range(n_rows):
        c=0
        a,c=partial_Pivot(a,i)
        count+=c
        pivot=a[i][i]
        for k in range(n_rows):
            if k!=i:
                sub=a[k][i]
                for l in range(n_rows):
                    a[k][l]-=(sub/pivot)*a[i][l]
    det=(-1)**count#To account for the change in sign of the determinant when the rows are swapped
    for i in range(n_rows):
        det*=a[i][i]
    return det      


#Finds the inverse of the matrix using crout's decomposition
def inverse_Matrix(a):
    n_rows=len(a)
    identity=[]
    for i in range(n_rows):
        row = [1 if (i==j) else 0 for j in range(n_rows)]
        identity.append(row)
    
    return(forwardbackwardsub_crouts(lu_crouts(a,identity)))

#Functions to find the number of rows and columns in a matrix
row_no = lambda mat:len(mat)
col_no = lambda mat:len(mat[0])

#Function to multiply two matrices
def multiply_Matrix(mat1,mat2):
    if col_no(mat1) == row_no(mat2):
        return [[sum(mat1[i][k]*mat2[k][j] for k in range(col_no(mat1)))
                 for j in range(col_no(mat2))] for i in range(row_no(mat1))]
    else:
        print("Indexes do not match")


#Function to find the LU decomposition of a matrix using doolittle's meathod
def lu_doolittle(a,b):
    ab=aug_Matrix(a,b)
    row=row_no(a)
    for i in range(row):
        ab,_=partial_Pivot(ab,i)
        for j in range(row):
            if i<=j:#Finding elements of upper triangular matrix
                f = sum(ab[i][k]*ab[k][j] for k in range(i))
                ab[i][j]=ab[i][j]-f
            else:#Finding elements of lower triangular matrix
                f = sum(ab[i][k]*ab[k][j] for k in range(j))
                ab[i][j]=(ab[i][j]-f)/ab[j][j]
    return(ab)

#Function to find the LU decomposition of a matrix using crout's meathod
def lu_crouts(a,b):
    ab=aug_Matrix(a,b)
    row=row_no(a)
    for i in range(row):
        ab,_=partial_Pivot(ab,i)
        for j in range(row):
            if j<=i:#Finding elements of lower triangular matrix
                f = sum(ab[i][k]*ab[k][j] for k in range(j))
                ab[i][j]=ab[i][j]-f
            else:#Finding elements of upper triangular matrix
                f = sum(ab[i][k]*a[k][j] for k in range(i))
                ab[i][j]=(ab[i][j]-f)/ab[i][i]
    return(ab)


#Finds the determinant of upper triangular matrix
def determinant_upper(a):
    row=row_no(a)
    det=1
    for i in range(row):
        det*=a[i][i]
    return det

#Does the forward and backward substitution to solve the equation or find inverse in case of doolittle's decomposition
def forwardbackwardsub_doolittle(ab):
    arow=row_no(ab)
    bcol=col_no(ab)-arow
    y = [[0 for y in range(bcol)] for x in range(arow)]
    for i in range(arow):#forward substitution
        for j in range(bcol):
            f = sum(ab[i][k]*y[k][j] for k in range(i))
            y[i][j]=(ab[i][j+arow]-f)
    x = [[0 for y in range(bcol)] for x in range(arow)]
    for i in range(arow-1, -1, -1):#backward substituion
        for j in range(bcol):
            f = sum(ab[i][k]*x[k][j] for k in range(i + 1, arow))
            x[i][j]=(y[i][j]-f)/ab[i][i]
    return x


#Does the forward and backward substitution to solve the equation or find inverse in case of crout's decomposition
def forwardbackwardsub_crouts(ab):
    arow=row_no(ab)
    bcol=col_no(ab)-arow
    y = [[0 for y in range(bcol)] for x in range(arow)]
    for i in range(arow):#forward substitution
        for j in range(bcol):
            f = sum(ab[i][k]*y[k][j] for k in range(i))
            y[i][j]=(ab[i][j+arow]-f)/ab[i][i]
    x = [[0 for y in range(bcol)] for x in range(arow)]
    for i in range(arow-1, -1, -1):#backward substitution
        for j in range(bcol):
            f = sum(ab[i][k]*x[k][j] for k in range(i+1, arow))
            x[i][j]=(y[i][j]-f)
    return x

#Function to find the LU decomposition of a matrix using cholesky's meathod
'''
Works only for Hermitian and positive definite matrices
'''
def lu_cholesky(a,b):
    ab=aug_Matrix(a,b)
    arow=row_no(ab)
    for i in range(arow):
        ab,_ = partial_Pivot(ab,i)
        for j in range(i+1):
            if i == j:
                f = sum(ab[j][k]**2 for k in range(j))
                ab[j][j] = float(sqrt(ab[j][j] - f))
            if i > j:
                f = sum(ab[i][k]*ab[j][k] for k in range(j))
                if ab[j][j] > 0 :
                    ab[i][j] = float((ab[i][j] - f)/ab[j][j])
                    ab[j][i]=ab[i][j]
    return ab

def forwardbackwardsub_cholesky(ab):
    arow=row_no(ab)
    bcol=col_no(ab)-arow
    y = [[0 for y in range(bcol)] for x in range(arow)]
    for i in range(arow):
        for j in range(bcol):
            f = sum(ab[i][k]*y[k][j] for k in range(i))
            y[i][j]=(ab[i][j+arow]-f)/ab[i][i]
    x = [[0 for y in range(bcol)] for x in range(arow)]
    for i in range(arow-1, -1, -1):
        for j in range(bcol):
            f = sum(ab[i][k]*x[k][j] for k in range(i+1, arow))
            x[i][j]=(y[i][j]-f)/ab[i][i]
    return x