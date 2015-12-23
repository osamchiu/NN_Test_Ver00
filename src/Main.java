public class Main {
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		double train_u[][] = new double[2][800];
		double train_y[][] = new double[1][800];
		double test_u[][] = new double[2][200];
		double test_y[][] = new double[1][200];
		
		Matrix matrix = new Matrix();
		NN nn = new NN(2, 1, 3);
		nn.Initialization();
		RandomCreateUserDefineMatrix(train_u,train_y);
		RandomCreateUserDefineMatrix(test_u,test_y);

		//matrix.ShowMatrix(train_u);
		//matrix.ShowMatrix(train_y);
		
		matrix.NormalizationMatrix(train_u);
		matrix.NormalizationMatrix(train_y);
		matrix.NormalizationMatrix(test_u);
		matrix.NormalizationMatrix(test_y);

		nn.TrainingBox(train_u, train_y);
		nn.TestingBox(test_u, test_y);
		
		//System.out.println("error: " + nn.error[0][0]);
		//nn.TrainingBox(test_u, test_y);
		
		//System.out.println("Out error: " + nn.error[0][0]);
		
		/*
		double[][] m1 = {	{0.2},
							{0.3}
						};
        double[][] m2 = {
        					{0.5}
        				};
        
		nn.forward(m1);
		nn.backward(m1, m2);
		double[][] m3 = {	{0.1361654},
							{0.521}
						};
		double[][] m4 = {
							{0.95874321}
						};

		nn.forward(m1);
		nn.backward(m1, m2);
		
		
		//nn.printNN();
		
		double[][] m1 = { { 5, 2 }, { 4, 2 }, { 52, 1 } };
        double[][] m2 = { { 1, 2 }, { 3, 4 }, { 5, 6 } };
        double[][] m3;
        Matrix matrix = new Matrix();
        m3 = matrix.SubtractionMatrix(m1, m2);
		matrix.ShowMatrix(m3);
		
		double[][] m1 = { { 5, 2 }, { 4, 2 }, { 52, 1 } };
        double[][] m2 = { { 1, 2 }, { 3, 4 }, { 5, 6 } };
        double[][] m3;
        //double[][] x = new double(m3);
        Matrix matrix = new Matrix();
        m3 = matrix.MultiplyingMatrix(m1, m2);
		matrix.ShowMatrix(m3);
		
		
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
				u[i][j] = Math.random()*2;
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
	
	double e = 0.001;
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
		RandomCreatParameterMatrixZeroToOne(this.sigma);
		RandomCreateParameterMatrix(this.w3);
	}
	public void RandomCreateParameterMatrix(double x[][]){
		for(int i = 0; i<x.length; i++){
			for(int j = 0; j<x[i].length; j++){
				x[i][j] = Math.random()*2-1;
			}
		}
	}
	public void RandomCreatParameterMatrixZeroToOne(double x[][]){
		for(int i = 0; i<x.length; i++){
			for(int j = 0; j<x[i].length; j++){
				x[i][j] = Math.random();
			}
		}
	}
	
	public void forward(double input[][]){
		Matrix matrix = new Matrix();
		//double[][] Qr = new double[this.numberOfOutput][1];
		
		this.fq = matrix.MultiplyingMatrix(this.w1, input);
		for(int i=0; i<this.Oq.length; i++){
			this.Oq[i][0] = (Math.exp(this.fq[i][0])-Math.exp(-1*this.fq[i][0]))/(Math.exp(this.fq[i][0])+Math.exp(-1*this.fq[i][0]));
		}

		this.Oq = matrix.SubtractionMatrix(this.fq, this.m);
		for(int i=0; i<this.Oq.length; i++){
			this.Oq[i][0] *= (-1 * this.Oq[i][0]); 
			this.Oq[i][0] /= Math.pow(this.sigma[i][0],2);
			this.Oq[i][0] = Math.exp(this.Oq[i][0]);
		}

		this.Or = matrix.MultiplyingMatrix(this.w3, this.Oq);

		//return null;
	}
	
	public void backward(double input[][], double output[][]){
		Matrix matrix = new Matrix();
		//System.out.println("Before error: " + this.error[0][0]);
		this.error = matrix.SubtractionMatrix(output, this.Or);
		//System.out.println("error: " + this.error[0][0]);
		//double[][] tempw1 = matrix.CopyMatrix(this.w1);
		double[][] tempw3 = matrix.CopyMatrix(this.w3);
		double[][] tempm = matrix.CopyMatrix(this.m);
		double[][] tempsigma = matrix.CopyMatrix(this.sigma);
		for(int i=0; i<this.numberOfW; i++){
			this.w3[0][i] += this.e * error[0][0] * this.Oq[i][0]; 
			this.m[i][0] += (-2 * this.e * this.error[0][0] * tempw3[0][i] * (this.fq[i][0] - tempm[i][0]) * this.Oq[i][0] / Math.pow(tempsigma[i][0],2));
			this.sigma[i][0] += (-2 * this.e * this.error[0][0] * tempw3[0][i] * Math.pow((this.fq[i][0]-tempm[i][0]), 2) * this.Oq[i][0] / Math.pow(tempsigma[i][0],3));;

			for(int j=0; j<this.w1[i].length; j++){
				this.w1[i][j] = (2 * this.e * this.error[0][0] * tempw3[0][i] * (this.fq[i][0] - tempm[i][0]) * this.Oq[i][0] * input[j][0]/ Math.pow(tempsigma[i][0],2));
			}

			/*
			for(int j=0; j<this.w1[i].length; j++){
				this.w1[i][j] = (-4 * this.error[0][0] * tempw3[0][i] * input[j][0] / Math.pow((Math.exp(this.fq[i][0]) + Math.exp(-1*this.fq[i][0])),2));
			}
			*/
		}
	}
		
	public void TrainingBox(double input[][], double output[][]){
		double[][] tempInput = new double[input.length][1];
		double[][] tempOutput = new double[output.length][1];
		double[][] tempError = new double[output.length][output[0].length];
		int epochs = 1000;
		
		for(int i=0; i<epochs; i++) {

			for (int j = 0; j < input[0].length; j++) {
				for (int k = 0; k < input.length; k++) {
					tempInput[k][0] = input[k][j];
				}
				for (int k = 0; k < output.length; k++) {
					tempOutput[k][0] = output[k][j];
				}
				forward(tempInput);
				backward(tempInput, tempOutput);
				//要修改 有爆破
				tempError[0][j] = this.error[0][0];
			}

			double sum = 0;
			for (int j = 0; j < tempError[0].length; j++) {
				sum += Math.pow(tempError[0][j], 2);
			}
			sum /= tempError[0].length;
			sum = Math.sqrt(sum);
			System.out.println("Training RMS(error): " + sum);
		}

		/*
		Matrix matrix = new Matrix();
		for(int j=0; j<input[0].length; j++){
			for(int k=0; k<input.length; k++){
				tempInput[k][0] = input[k][j];
			}
			for(int k=0; k<output.length; k++){
				tempOutput[k][0] = output[k][j];
			}
			//matrix.ShowMatrix(tempInput);
			//matrix.ShowMatrix(tempOutput);
			forward(tempInput);
			backward(tempInput, tempOutput);
		}
		*/
	}
	public void TestingBox(double input[][], double output[][]){

		double[][] tempInput = new double[input.length][1];
		double[][] tempOutput = new double[output.length][1];
		double[][] tempError = new double[output.length][output[0].length];
		double sum = 0;
		for (int j = 0; j < input[0].length; j++) {
			for (int k = 0; k < input.length; k++) {
				tempInput[k][0] = input[k][j];
			}
			forward(tempInput);
			tempError[0][j] = output[0][j] - this.Or[0][0];
			//System.out.println("tempError = " + tempError[0][j]);
			sum += Math.pow(tempError[0][j], 2);
			System.out.println("Test RMS(error): " + Math.sqrt(sum/j));

			//sum /= j;
			//sum = Math.sqrt(sum);
			//System.out.println("Test RMS(error): " + sum);

		}

	}
	public void printNN(){
		Matrix matrix = new Matrix();

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
