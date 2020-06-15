package moead;

import Utils.EDC;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;


public class Initializer {

	private static double timeLimit; //时间约束
	private static int FPGA_Square_Limit; //FPGA面积计算

	public Initializer(double timeLimit,int FPGA_Square_Limit){
		this.timeLimit = timeLimit;
		this.FPGA_Square_Limit = FPGA_Square_Limit;
	}

	public Initializer(){
		timeLimit = 500.0;
		FPGA_Square_Limit = 4500;
	}

	public static double getTimeLimit() {
		return timeLimit;
	}

	public void setTimeLimit(double timeLimit) {
		this.timeLimit = timeLimit;
	}

	public int getFPGA_Square_Limit() {
		return FPGA_Square_Limit;
	}

	public void setFPGA_Square_Limit(int FPGA_Square_Limit) {
		this.FPGA_Square_Limit = FPGA_Square_Limit;
	}



	//产生二维权重向量 N是产生的权重个数
	public static double[][] generateWeightVectors(int N) {
		// fixed objective size to 2 
		double[][] weightVectors = new double[N][2];//[N][0]:X轴 [N][1]:Y轴

		double slice = 1.0 / (N - 1);//使得向量的XY分量均在0-1之间
		for (int i = 0; i < N; i++) {
			weightVectors[i][0] = i * slice;//X轴分量
			weightVectors[i][1] = (1 - i * slice);//Y轴分量
		}

		return weightVectors;
	}

	//获得邻居
	public static int[][] getNeighbors(double[][] weightVectors, int neighborSize) {
		int N = weightVectors.length;//二维数组行数
		double[][] distances = getDistances(weightVectors);

		int[][] neighbors = new int[N][neighborSize];
		for (int i = 0; i < N; i++) {
			neighbors[i] = getSmallestValues(distances[i], neighborSize);
		}
		return neighbors;
	}

	//产生随机种群大小 genome:基因组
	public static double[][] getRandomPopulation(int N, int genomeSize) {
		Random rand = new Random();
		double[][] population = new double[N][genomeSize];
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < genomeSize; j++) {
				population[i][j] = rand.nextDouble();//返回一个大于或等于 0.0 且小于 1.0 的随机浮点数。
			}
		}
		return population;
	}

	//计算种群中每个点的关于两个目标函数的值
	public static double[][] computeFunctionValues(double[][] population) {
		double[][] values = new double[population.length][2];
		for (int i = 0; i < population.length; i++) {
			values[i][0] = Functions.f1(population[i]);//values[i][0]存储第一个目标函数值
			values[i][1] = Functions.f2(population[i]);//values[i][1]存储第二个目标函数值
		}
		return values;
	}

	//获取参考点
	public static double[] getReferencePoint(double[][] functionValues) {
		double[] referencePoint = new double[2];

		List<Double> f1Values = new ArrayList<Double>();
		List<Double> f2Values = new ArrayList<Double>();
		for (int i = 0; i < functionValues.length; i++) {
			f1Values.add(functionValues[i][0]);//存储所有点的第一个目标函数值
			f2Values.add(functionValues[i][1]);//存储所有点的第二个目标函数值
		}
		
		referencePoint[0] = Collections.min(f1Values);//找出所有点中第一个目标函数值最小值
		referencePoint[1] = Collections.min(f2Values);//找出所有点中第二个目标函数值最小值
		return referencePoint;//返回参考点
	}

	//获取距离
	private static double[][] getDistances(double[][] weightVectors) {
		int N = weightVectors.length;
		double[][] distances = new double[N][N];

		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				double distance = Math
						.sqrt(Math.pow(weightVectors[i][0]
								- weightVectors[j][0], 2)
								+ Math.pow(weightVectors[i][1]
										- weightVectors[j][1], 2));
				distances[i][j] = distance;
				distances[j][i] = distance;
			}
		}
		return distances;
	}
	
	private static int[] getSmallestValues(double[] distances, int neighborSize) {
		List<NeighborHelper> sortingHelper = new ArrayList<NeighborHelper>();
		Initializer initializer = new Initializer();
		for (int i = 0; i < distances.length; i++) {
			NeighborHelper n = initializer.new NeighborHelper(distances[i], i);
			sortingHelper.add(n);
		}
		Collections.sort(sortingHelper);//按distance从小到大的顺序排列，对应相应的index
		
		
		int[] neighbors = new int[neighborSize];
		for (int i = 0; i < neighborSize; i++) {
			neighbors[i] = sortingHelper.get(i).index;//把index存进去，distance越小的排的越前
		}
		return neighbors;
	}

	//实现Comparable接口才能使用Collection接口进行排序
	private class NeighborHelper implements Comparable<NeighborHelper>{
		double distance;
		int index;
		public NeighborHelper(double distance, int index) {
			this.distance = distance;
			this.index = index;
		}

		//目的是为了更改Collection.sort()的排序方式(返回值大于0表示正序(升序)，小于0表示逆序)
		//指定按distance正序排列
		@Override
		public int compareTo(NeighborHelper o) {
			double diff = this.distance - o.distance;
			return diff == 0 ? 0 : (diff < 0 ? -1 : 1);
		}
	}




	//产生EDC随机种群大小
	/*
	初始化的种群应该是一个Poplation[][] 二维数组，每个元素都是一个EDC实例。
	Poplation的每一行(Poplation[k])都是一组可能的EDC消息排列，每一行都有对应的一个AFC和HDFC值
	 */
	public static EDC[][] getRandomEDCPopulation(int N,int genomeSize){
		Random random = new Random();
		EDC[][] edcPop = new EDC[N][genomeSize];

		double total_Time;
		int total_FPGA_Square;

		int soh;
		for (int i =0;i<N;i++){
			for (int j=0;j<genomeSize;j++){

				soh = random.nextInt(2);
//				edcPop[i][j] = new EDC(random.nextDouble());
				edcPop[i][j] = generateEDCduetoTable(soh,random);

			}

			total_Time = Functions.Total_Time(edcPop[i]);
			total_FPGA_Square = Functions.Total_FPGA_SQUARE(edcPop[i]);

			//新产生的EDC组 是否满足约束条件 如果不满足约束条件则把该Row重新覆盖
			if (total_Time > timeLimit || total_FPGA_Square > FPGA_Square_Limit){
				i--;
			}

		}
		return edcPop;
	}



	//计算种群中每个点的关于两个目标函数的值
	public static double[][] computeEDCFunctionValues(EDC[][] edc){
		double[][] values = new double[edc.length][4];
		for(int i =0;i<edc.length;i++){
			values[i][0] = Functions.Average_Coverage(edc[i]);//values[i][0]存储AFC
			values[i][1] = Functions.Calculate_std(edc[i]);//values[i][1]存储异构度
			values[i][2] = Functions.Total_Time(edc[i]);//values[i][2]存储总时间
			values[i][3] = Functions.Total_FPGA_SQUARE(edc[i]);//values[i][3]存储总FPGA使用面积
		}
		return values;
	}


	public static EDC generateEDCduetoTable(int soh,Random random){
		EDC edc = null;
		int num;
		do {
			num = random.nextInt(11);//产生0-10的整数
		}while (num < 1);

		switch (soh){
			case 0:// 选择软件实现方式
				switch (num){
					case 1:
						edc = new EDC(1,0.8716,19.5,soh);
						break;
					case 2:
						edc = new EDC(2,0.8719,22,soh);
						break;
					case 3:
						edc = new EDC(3,0.9126,24.5,soh);
						break;
					case 4:
						edc = new EDC(4,0.924,26,soh);
						break;
					case 5:
						edc = new EDC(5,0.9249,28.5,soh);
						break;
					case 6:
						edc = new EDC(6,0.9291,30,soh);
						break;
					case 7:
						edc = new EDC(7,0.9337,32.5,soh);
						break;
					case 8:
						edc = new EDC(8,0.9457,34,soh);
						break;
					case 9:
						edc = new EDC(9,0.9732,35.5,soh);
						break;
					case 10:
						edc = new EDC(10,1.0,41,soh);
						break;
				}
			break;

			case 1:// 选择硬件

				switch (num){
					case 1:
						edc = new EDC(1,0.8716,1.95,soh,600);
						break;
					case 2:
						edc = new EDC(2,0.8719,2.2,soh,614);
						break;
					case 3:
						edc = new EDC(3,0.9126,2.45,soh,724);
						break;
					case 4:
						edc = new EDC(4,0.924,2.6,soh,609);
						break;
					case 5:
						edc = new EDC(5,0.9249,2.85,soh,611);
						break;
					case 6:
						edc = new EDC(6,0.9291,3.0,soh,627);
						break;
					case 7:
						edc = new EDC(7,0.9337,3.25,soh,727);
						break;
					case 8:
						edc = new EDC(8,0.9457,3.4,soh,827);
						break;
					case 9:
						edc = new EDC(9,0.9732,3.55,soh,925);
						break;
					case 10:
						edc = new EDC(10,1.0,4.1,soh,1168);
						break;
				}

			break;

			default:

		}


		return edc;

	}






}
