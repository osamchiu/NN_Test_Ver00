public class Main {
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		double train_u[][] = new double[2][80];
		double train_y[][] = new double[1][80];
		double test_u[][] = new double[2][20];
		double test_y[][] = new double[1][20];

		
		RandomCreateUserDefineMatrix(train_u,train_y);
		RandomCreateUserDefineMatrix(test_u,test_y);
		
		
		double[][] m1 = { { 1, 2, 3 }, { 4, 5, 6 } };
        double[][] m2 = { { 1, 2 }, { 3, 4 }, { 5, 6 } };
        double[][] m3;
        //double[][] x = new double(m3);
        Matrix matrix = new Matrix();
        m3 = matrix.MultiplyingMatrix(m1, m2);
		matrix.ShowMatrix(m3);
		
		/*
		NormalizationMatrix(u[0]);
		NormalizationMatrix(u[1]);
		NormalizationMatrix(y[0]);

		RandomCreateParameterMatrix(w1);
		RandomCreateParameterMatrix(m);
		RandomCreateParameterMatrix(sigma);
		RandomCreateParameterMatrix(w3);
		
		for(int i=0; i<fq.length; i++){
			for(int j=0; j<fq[i].length; j++){
				fq[i][j] = 0;
			}
		}
		
		*/
		
		
		
		
		/*
		System.out.println("u");
		ShowMatrix(u);
		System.out.println("---------------------------------------------------------------------------------------------------------------------------------------------------");
		System.out.println("y");
		ShowMatrix(y);
		System.out.println("---------------------------------------------------------------------------------------------------------------------------------------------------");
		System.out.println("w1");
		ShowMatrix(w1);
		System.out.println("---------------------------------------------------------------------------------------------------------------------------------------------------");
		System.out.println("m");
		ShowMatrix(m);
		System.out.println("---------------------------------------------------------------------------------------------------------------------------------------------------");
		System.out.println("sigma");
		ShowMatrix(sigma);
		System.out.println("---------------------------------------------------------------------------------------------------------------------------------------------------");
		System.out.println("w3");
		ShowMatrix(w3);
		System.out.println("---------------------------------------------------------------------------------------------------------------------------------------------------");
		
		System.out.println("TransferMatrix u");
		double test[][] = TransferMatrix(u);
		ShowMatrix(test);
		System.out.println("---------------------------------------------------------------------------------------------------------------------------------------------------");
		*/
		
	}
	
	public static void RandomCreateUserDefineMatrix(double u[][], double y[][]){
		for(int i = 0; i<y.length; i++){
			for(int j = 0; j<y[i].length; j++){
				y[i][j] = 0;		
			}
		}
		
		for(int i = 0; i<u.length; i++){
			for(int j = 0; j<u[i].length; j++){
				u[i][j] = Math.random()*10;
				y[0][j] += u[i][j]*u[i][j];
			}
		}

	}
	
	

	
}

class Matrix{
	public void ShowMatrix(double x[][]){
		for(int i = 0; i < x.length; i++){
			for(int j = 0;j<x[i].length; j++){
				System.out.print(x[i][j] + "	");
			}
			System.out.println();
		}
	}
	
	public double[][] TransferMatrix(double x[][]){
		double ret[][] = new double[x[0].length][x.length];
		for(int i=0; i<x.length; i++){
			for(int j=0; j<x[i].length; j++){
				ret[j][i] = x[i][j];
			}
		}
		return ret;
	}
	
	public void NormalizationMatrix(double x[][]){
		double max = 0, min = 0;
		for(int i=0; i<x.length; i++){
			for(int j=0; j<x[i].length; j++){
				if(x[i][j] > max){max = x[i][j];}
				if(x[i][j] < min){min = x[i][j];}
			}
			
			for(int j=0; j<x[i].length; j++){
			x[i][j] = ((x[i][j] - min)/(max - min)*2)-1;
			}
		}
	}
	
	public double[][] MultiplyingMatrix(double x[][], double y[][]){
		double[][] ret = new double[x.length][y[0].length];
		double sum = 0;
        for (int i = 0;i< ret.length; i++) {
            for (int j = 0; j < ret[i].length; j++) {
                for (int s = 0; s < y.length; s++) {
                    sum += x[i][s] * y[s][j];
                }
                ret[i][j] = sum;
                sum = 0;
            }
        }
        return ret;
	}
	public double[][] SubtractionMatrix(double x[][], double y[][]){
		if((x.length == y.length)&&(x[0].length == y[0].length)){
			double[][] ret = new double[x.length][x[0].length];
			for(int i=0; i<ret.length; i++){
				for(int j=0; j<ret[i].length; j++){
					ret[i][j] = x[i][j] - y[i][j];
				}
			}
			return ret;
		}
		else{
			System.out.println("SubtractionMatrix matrix size doesn't equal!!!");
			return null;
			
		}
	}
	public double[][] CopyMatrix(double x[][]){
		double[][] ret = new double[x.length][x[0].length];
		for(int i=0; i<x.length; i++){
			for(int j=0; j<x[i].length; j++){
				ret[i][j] = x[i][j];
			}
		}
		return ret;	
	}
}

class NN{
	int numberOfW;
	int numberOfInput;
	int numberOfOutput;
	
	double w1[][];
	double m[][];
	double sigma[][];
	double w3[][];
	
	double fq[][];
	double Oq[][];
	double Or[][];
	double error[][];
	
	double e = 0.1;
	NN(int numberOfInput, int numberOfOutput, int numberOfW){
		this.numberOfW = numberOfW;
		this.numberOfInput = numberOfInput;
		this.numberOfOutput = numberOfOutput;

		this.w1 = new double[this.numberOfW][this.numberOfInput];
		this.m = new double[this.numberOfW][1];
		this.sigma = new double[this.numberOfW][1];
		this.w3 = new double[this.numberOfOutput][this.numberOfW];
		
		this.fq = new double[this.numberOfW][1];
		this.Oq = new double[this.numberOfW][1];
		this.Or = new double[this.numberOfOutput][1];
		this.error = new double[this.numberOfOutput][1];
	}
	
	
	public void Initialization(){
		RandomCreateParameterMatrix(this.w1);
		RandomCreateParameterMatrix(this.m);
		RandomCreateParameterMatrix(this.sigma);
		RandomCreateParameterMatrix(this.w3);
	}
	public void RandomCreateParameterMatrix(double x[][]){
		for(int i = 0; i<x.length; i++){
			for(int j = 0; j<x[i].length; j++){
				x[i][j] = Math.random()*2-1;
			}
		}
	}
	
	public double[][] forward(double input[][]){
		Matrix matrix = new Matrix();
		double[][] Qr = new double[this.numberOfOutput][1];
		
		this.fq = matrix.MultiplyingMatrix(this.w1, input);
		this.Oq = matrix.SubtractionMatrix(this.fq, this.m);
		for(int i=0; i<this.Oq.length; i++){
			this.Oq[i][0] *= (-1 * this.Oq[i][0]); 
			this.Oq[i][0] /= (this.sigma[i][0] * this.sigma[i][0]);
			this.Oq[i][0] = Math.exp(this.Oq[i][0]);
		}
		this.Or = matrix.MultiplyingMatrix(this.w3, this.Oq);
		return null;
	}
	
	public void backward(double input[][], double output[][]){
		Matrix matrix = new Matrix();
		this.error = matrix.SubtractionMatrix(output, this.Or);
		double[][] tempw1 = matrix.CopyMatrix(this.w1);
		double[][] tempw3 = matrix.CopyMatrix(this.w3);
		double[][] tempm = matrix.CopyMatrix(this.m);
		double[][] tempsigma = matrix.CopyMatrix(this.sigma);
		for(int i=0; i<this.numberOfW; i++){
			this.w3[i][0] += this.e * error[0][0] * this.Oq[i][0]; 
			this.m[i][0] += (-2 * this.e * this.error[0][0] * tempw3[i][0] * (this.fq[i][0] - this.m[i][0]) * this.Oq[i][0] / Math.pow(this.sigma[i][0],2));
			this.sigma[i][0] += (-2 * this.e * this.error[0][0] * tempw3[i][0] * Math.pow((this.fq[i][0]-tempm[i][0]), 2) * this.Oq[i][0] / Math.pow(this.sigma[i][0],3));;
			this.w1[i][0] = (2 * this.e * this.error[0][0] * tempw3[i][0] * (this.fq[i][0] - this.m[i][0]) * this.Oq[i][0] / Math.pow(this.sigma[i][0],2));
		}
	}
		
	public void TrainingBox(double input[][], double output[][]){
		double[][] temp = new double[input.length][1];
		
	}
	
	public void printNN(){
		Matrix matrix = new Matrix();
		/*
		System.out.println("input");
		matrix.ShowMatrix(this.input);
		System.out.println("---------------------------------------------------------------------------------------------------------------------------------------------------");
		
		System.out.println("output");
		ShowMatrix(this.output);
		System.out.println("---------------------------------------------------------------------------------------------------------------------------------------------------");
		*/
		System.out.println("w1");
		matrix.ShowMatrix(this.w1);
		System.out.println("---------------------------------------------------------------------------------------------------------------------------------------------------");
		
		System.out.println("m");
		matrix.ShowMatrix(this.m);
		System.out.println("---------------------------------------------------------------------------------------------------------------------------------------------------");
		
		System.out.println("sigma");
		matrix.ShowMatrix(this.sigma);
		System.out.println("---------------------------------------------------------------------------------------------------------------------------------------------------");
		
		System.out.println("w3");
		matrix.ShowMatrix(this.w3);
		System.out.println("---------------------------------------------------------------------------------------------------------------------------------------------------");
		
		System.out.println("fq");
		matrix.ShowMatrix(this.fq);
		System.out.println("---------------------------------------------------------------------------------------------------------------------------------------------------");
		
		System.out.println("Oq");
		matrix.ShowMatrix(this.Oq);
		System.out.println("---------------------------------------------------------------------------------------------------------------------------------------------------");
		
		System.out.println("Or");
		matrix.ShowMatrix(this.Or);
		System.out.println("---------------------------------------------------------------------------------------------------------------------------------------------------");
		
		System.out.println("error");
		matrix.ShowMatrix(this.error);
		System.out.println("---------------------------------------------------------------------------------------------------------------------------------------------------");
		
	}
}
