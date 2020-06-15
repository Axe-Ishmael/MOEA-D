package moead;

import Utils.EDC;

import java.util.List;

public class Functions {


	public static double f1(List<Double> x) {
		double result = x.get(0);
		int n = x.size();
		
		double sum = 0;
		double count = 0; // length of J
		for (int j = 3; j < n; j += 2) {
			double pow = 0.5 * (1.0 + (3 * (j - 2)) / (n + 2));
			double term = Math.pow(x.get(j) - Math.pow(x.get(0), pow), 2);
			sum += term;
			count++;
		}
		
		return result + (2 / count) * sum;
	}

	public static double f2(List<Double> x) {
		double result = 1 - Math.sqrt(x.get(0));
		int n = x.size();
		
		double sum = 0;
		double count = 0; // length of J
		for (int j = 2; j < n; j += 2) {
			double pow = 0.5 * (1.0 + (3 * (j - 2)) / (n + 2));
			double term = Math.pow(x.get(j) - Math.pow(x.get(0), pow), 2);
			sum += term;
			count++;
		}
		
		return result + (2 / count) * sum;
	}
	
	public static double f1(double[] x) {
		double result = x[0];
		int n = x.length;
		
		double sum = 0;
		double count = 0; // length of J
		for (int j = 3; j < n; j += 2) {
			double pow = 0.5 * (1.0 + (3 * (j - 2)) / (n + 2));
			double term = Math.pow(x[j] - Math.pow(x[0], pow), 2);
			sum += term;
			count++;
		}
		
		return result + (2 / count) * sum;
	}

	public static double f2(double[] x) {
		double result = 1 - Math.sqrt(x[0]);
		int n = x.length;
		
		double sum = 0;
		double count = 0; // length of J
		for (int j = 2; j < n; j += 2) {
			double pow = 0.5 * (1.0 + (3 * (j - 2)) / (n + 2));
			double term = Math.pow(x[j] - Math.pow(x[0], pow), 2);
			sum += term;
			count++;
		}
		
		return result + (2 / count) * sum;
	}
	
	public static double f1Norm(double[] x, double min, double max) {
		return f1(x) / (max - min);
	}

	public static double f2Norm(double[] x, double min, double max) {
		return f2(x) / (max - min);		
	}

//---------------------------------------------EDC Part-------------------------------------------------


	public static double MFT(EDC[] edc){
		double sum = 0.0;
		double FTO = 0.99;
		for(int i = 0;i<edc.length;i++){
			sum += (1-edc[i].getError_Coverage())/(1-FTO);
		}
		return sum/edc.length;
	}

	public static double HDFT(EDC[] edc){
		double FTO = 0.99;  //目标错误覆盖率
		double average = MFT(edc);
		double sum = 0;
		double FT = 0;
		for (int i = 0;i<edc.length;i++){
			FT = (1 - edc[i].getError_Coverage()) / (1 - FTO);
			sum += (FT - average)*(FT - average);
		}
		sum = Math.sqrt(sum / edc.length);
		return sum;
	}

	public static double obj(EDC[] edc){
		double sum = 0.9*MFT(edc) + 0.1*HDFT(edc);
		return sum;
	}

	//求解平均错误覆盖率
	public static double Average_Coverage( EDC[] edc)
	{
		double sum = 0.0;
		for (int i = 0; i < edc.length; i++)
		{
			sum += edc[i].getError_Coverage();
		}
		return sum / (edc.length);
	}

	//求解异构度
	public static double Calculate_std(EDC[] edc)
	{
		double average = Average_Coverage(edc);
		double sum = 0;
		for (int i = 0; i < edc.length; i++)
		{
			sum += (edc[i].getError_Coverage() - average)*(edc[i].getError_Coverage() - average);
		}
		sum = Math.sqrt(sum / edc.length);
		return sum;
	}

	public static double Total_Time(EDC[] edc){

		double total_time = 0;
		for (int i = 0 ;i<edc.length;i++){
			total_time += edc[i].getEDC_Time();
		}

		return total_time;
	}


	public static int Total_FPGA_SQUARE(EDC[] edc){

		int total_FPGA_Square = 0;
		for (int i = 0;i<edc.length;i++){
			total_FPGA_Square += edc[i].getFPGA_Sqare();
		}

		return total_FPGA_Square;

	}


}
