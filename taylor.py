from sage.all import RIF, vector, zero_vector, sin, cos, exp

class Taylor:
    def __init__(self, val , degree, variable = False, field = RIF):
        
        if isinstance(val, Taylor):
            self.val = val.val
            self.degree = val.degree
        elif isinstance(val, list):
            self.val = vector(field, val)
            self.degree = degree
        else:
            self.degree = degree
            self.val = zero_vector(field, degree+1)
            self.val[0] = field( val )
            if variable:
                self.val[1]=field(1.0)
        
    def __add__(self, other):
        if isinstance(other, Taylor) is False:
            other = Taylor(other, degree = self.degree)
        
        if self.degree == other.degree:
            return Taylor((self.val+other.val).list(), self.degree)
        else:
            raise RuntimeError('The two series have different degree')

    def __radd__(self, other):
        if isinstance(other, Taylor) is False:
            other = Taylor(other,degree = self.degree)
        
        if self.degree != other.degree:
            raise RuntimeError('The two series have different degree')
        else:
            return Taylor((self.val+other.val).list(), self.degree)
    
    def __sub__(self, other):
        if isinstance(other, Taylor) is False:
            other = Taylor(other, degree = self.degree)
        
        if self.degree == other.degree:
            return Taylor((self.val-other.val).list(), self.degree)
        else:
            raise RuntimeError('The two series have different degree')

    def __rsub__(self, other):
        if isinstance(other, Taylor) is False:
            other = Taylor(other,degree = self.degree)
        
        if self.degree != other.degree:
            raise RuntimeError('The two series have different degree')
        else:
            return Taylor((other.val-self.val).list(), self.degree)
    
    
    
    def __mul__(self, other):
        if isinstance(other, Taylor) is False:
            other = Taylor(other,degree = self.degree)
        
        if self.degree != other.degree:
            raise RuntimeError('The two series have different degree')
        else:
            val = zero_vector(self.val[0].parent(), self.degree+1)
            for k in range(0,self.degree+1):
                for i in range(k+1):
                    val[k] += self.val[i]*other.val[k-i]
            return Taylor(val.list(),self.degree)
        
    def __rmul__(self, other):
        if isinstance(other, Taylor) is False:
            other = Taylor(other,degree = self.degree)
        if self.degree != other.degree:
            raise RuntimeError('The two series have different degree')
        else:
            val = zero_vector(self.val[0].parent(), self.degree+1)
            for k in range(0,self.degree+1):
                for i in range(k+1):
                    val[k] += self.val[i]*other.val[k-i]
            return Taylor(val.list(),self.degree)
        
    def exp(self):
        val = zero_vector(self.val[0].parent(), self.degree+1)
        val[0] = exp(self.val[0])
        for k in range(1,self.degree+1):
            for i in range(1, k+1):
                val[k] += i*self.val[i]*val[k-i]
            val[k]=val[k]/k
        return Taylor(val.list(),self.degree)
   
    def sin(self):
        val_sin = zero_vector(self.val[0].parent(), self.degree+1)
        val_cos = zero_vector(self.val[0].parent(), self.degree+1)
        
        val_sin[0] = sin(self.val[0])
        val_cos[0] = cos(self.val[0])
        
        for k in range(1,self.degree+1):
            for i in range(1, k+1):
                val_sin[k] += i*self.val[i]*val_cos[k-i]
                val_cos[k] += i*self.val[i]*val_sin[k-i]
            val_sin[k]=val_sin[k]/k
            val_cos[k]=-1*val_cos[k]/k
        return Taylor(val_sin.list(),self.degree)
    
    def cos(self):
        val_sin = zero_vector(self.val[0].parent(), self.degree+1)
        val_cos = zero_vector(self.val[0].parent(), self.degree+1)
        
        val_sin[0] = sin(self.val[0])
        val_cos[0] = cos(self.val[0])
        
        for k in range(1,self.degree+1):
            for i in range(1, k+1):
                val_sin[k] += i*self.val[i]*val_cos[k-i]
                val_cos[k] += i*self.val[i]*val_sin[k-i]
            val_sin[k]=val_sin[k]/k
            val_cos[k]=-1*val_cos[k]/k
        return Taylor(val_cos.list(),self.degree)
