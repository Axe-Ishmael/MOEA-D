package moead;

import Utils.EDC;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

public class Moead {

	private double minF1 = Double.MAX_VALUE;
	private double maxF1 = -Double.MAX_VALUE;
	private double minF2 = Double.MAX_VALUE;
	private double maxF2 = -Double.MAX_VALUE;

	public List<List<Double>> moead(int iterations, int N, int neighborSize, int genomeSize) {
		double[][] weightVectors = Initializer.generateWeightVectors(N);//产生N个权重向量
		
		// 1.1
		// Dosen't work as expected, so was commented out
		// Instead full list is returned
		//List<List<Double>> EP = new ArrayList<List<Double>>(); // external population
		// 1.2
		int[][] Neighbors = Initializer.getNeighbors(weightVectors, neighborSize); // neighbors 获得邻居
		// 1.3
//		double[][] population = Initializer.getRandomPopulation(N, genomeSize);
		EDC[][] population = Initializer.getRandomEDCPopulation(N,genomeSize);
//		double[][] functionValues = Initializer.computeFunctionValues(population);
		double[][] functionValues = Initializer.computeEDCFunctionValues(population);
		// 1.4
		double[] refPoint = Initializer.getReferencePoint(functionValues); // reference point

		int count = 0;
		boolean end = false;
		while (!end) {

			updateMinMax(functionValues);
			for (int i = 0; i < N; i++) {
				// 2.1
				EDC[] newSolution = reproduce(population, Neighbors[i], i);
				// 2.2 - no repair needed
				// 2.3
				updateReferencePoint(refPoint, newSolution, weightVectors[i]);
				// 2.4
				updateNeighborhood(population, functionValues, Neighbors[i], weightVectors, refPoint, newSolution);
				// 2.5
//				Dosen't work as expected, so was commented out
//				updateEP(EP, y, weightVectors[i], z);
			}
			
			if (count % 100 == 0) {
				System.out.println(count);
			}
			if (++count >= iterations) {
				end = true;
			}
		}

		List<List<Double>> result = new ArrayList<List<Double>>();
		for (int i = 0; i < functionValues.length; i++) {//functionValues.length = Population.length(EDC组数)
			List<Double> l = new ArrayList<Double>();
			double val1 = (functionValues[i][0] - refPoint[0]) / (maxF1 - minF1);
			double val2 = (functionValues[i][1] - refPoint[1]) / (maxF2 - minF2);
			l.add(val1);
			l.add(val2);
			result.add(l);
			System.out.println(functionValues[i][0] + " " + functionValues[i][1]);
		}
		
		return result;
		
//		for (List<Double> l : EP) {
//			double val1 = l.get(0);
//			double val2 = l.get(1);
//			double weightVector1 = l.get(2);
//			double weightVector2 = l.get(3);
//			
//			val1 = (val1 - z[0]) / (maxF1 - minF1);
//			val2 = (val2 - z[1]) / (maxF2 - minF2);
//			
//			l.set(0, val1);
//			l.set(1, val2);
//			l.remove(3);
//			l.remove(2);
//		}
//		return EP;
	}

	//找出functionValues中两个目标函数的最大值和最小值
	private void updateMinMax(double[][] functionValues) {
		List<Double> f1Values = new ArrayList<Double>();
		List<Double> f2Values = new ArrayList<Double>();
		for (int i = 0; i < functionValues.length; i++) {
			f1Values.add(functionValues[i][0]);
			f2Values.add(functionValues[i][1]);
		}
		
		minF1 = Collections.min(f1Values);
		maxF1 = Collections.max(f1Values);
		minF2 = Collections.min(f2Values);
		maxF2 = Collections.max(f2Values);
	}

	//(?)变异操作
//	private double[] reproduce(double[][] population, int[] neighbors, int index) {
//		Random rand = new Random();
//		int parent1Index = neighbors[rand.nextInt(neighbors.length)];//生成一个范围在0~neighbors.length（不包含neighbors.length）内的任意正整数
//		int parent2Index = neighbors[rand.nextInt(neighbors.length)];
//		double[] parent1 = population[parent1Index];//从population中选择一行
//		double[] parent2 = population[parent2Index];
//
//		double[] newSolution = new double[parent1.length];
//		for (int i = 0; i < parent1.length; i++) {
//			if (rand.nextDouble() < 0.01) {
//				// mutate (see paper p. 718)  变异操作
//				newSolution[rand.nextInt(parent1.length)] = rand.nextDouble();
//			}
//			else {
//				newSolution[i] = (parent1[i] + parent2[i]) / 2;
//			}
//		}
//
//		return newSolution;
//	}

	private EDC[] reproduce(EDC[][] population, int[] neighbors, int index) {
		Random rand = new Random();
		int parent1Index = neighbors[rand.nextInt(neighbors.length)];
		int parent2Index = neighbors[rand.nextInt(neighbors.length)];
		EDC[] parent1 = population[parent1Index];
		EDC[] parent2 = population[parent2Index];

		EDC[] newSolution = new EDC[parent1.length];
		for (int i = 0; i < parent1.length; i++) {
			newSolution[i] = new EDC();
			if (rand.nextDouble() < 0.01) {
				// mutate (see paper p. 718)
				newSolution[rand.nextInt(parent1.length)] = new EDC(rand.nextDouble());
			}
			else {
				newSolution[i] = new EDC((parent1[i].getError_Coverage() + parent2[i].getError_Coverage()) / 2);
			}
		}

		return newSolution;
	}

	//更新参考点
	private void updateReferencePoint(double[] refPoint, EDC[] newSolution, double[] weightVector) {
		double f1Value = Functions.Average_Coverage(newSolution);
		double f2Value = Functions.Calculate_std(newSolution);
		refPoint[0] = Math.min(refPoint[0], f1Value);
		refPoint[1] = Math.min(refPoint[1], f2Value);
	}

	//更新邻居
	private void updateNeighborhood(EDC[][] population,
			double[][] functionValues, int[] neighbors,
			double[][] weightVectors, double[] refPoint, EDC[] newSolution) {

		double y1Val = Functions.Average_Coverage(newSolution);
		double y2Val = Functions.Calculate_std(newSolution);

		for (int i = 0; i < neighbors.length; i++) {
			int index = neighbors[i];
//
			double gNeighbor = computeMaxCombinedValues(functionValues[index][0], functionValues[index][1], weightVectors[index], refPoint);
			double gY = computeMaxCombinedValues(y1Val, y2Val, weightVectors[index], refPoint);
			if (gY <= gNeighbor) {
				population[index] = newSolution;
				functionValues[index][0] = y1Val;
				functionValues[index][1] = y2Val;
			}
		}
	}

	//需要改
	private double computeMaxCombinedValues(double f1Val, double f2Val, double[] weightVector,
			double[] refPoint) {
//		First possibility: use max
		return Math.max(
				(f1Val - refPoint[0]) / (maxF1 - minF1) * weightVector[0],
				(f2Val - refPoint[1]) / (maxF2 - minF2) * weightVector[1]);

		//Second possibility: don't use max but instead addition
//		return ((f1Val - refPoint[0]) / (maxF1 - minF1) * weightVector[0]) + ((f2Val - refPoint[1]) / (maxF2 - minF2) * weightVector[1]);

		//-------------------------------------------------------------------------------
//		return Math.max(
//				(f1Val - refPoint[0]) / maxF1  * weightVector[0],
//				(f2Val - refPoint[1]) / maxF2  * weightVector[1]);
	}

//	Dosen't work as expected, so was commented out
//	private void updateEP(List<List<Double>> EP, double[] y, double[] weightVector, double[] z) {
//		// remove all vectors dominated by F(y)
//		// check if any vector in EP dominates F(y)
//		double f1Val = Functions.f1(y);
//		double f2Val = Functions.f2(y);
//		boolean yIsDominated = false;
//		for (Iterator<List<Double>> it = EP.iterator(); it.hasNext();) {
//			List<Double> epVector = it.next();
//
//			double[] yValue = new double[2];
//			yValue[0] = f1Val;
//			yValue[1] = f2Val;
//
//			double[] epValue = new double[2];
//			epValue[0] = epVector.get(0);
//			epValue[1] = epVector.get(1);
//
//			if (dominates(yValue, epValue)) {
//				it.remove();
//			}
//			if (dominates(epValue, yValue)) {
//				yIsDominated = true;
//			}
//			if (yValue[0] == epValue[0] && yValue[1] == epValue[1]) {
//				yIsDominated = true;
//			}
//		}
//		if (!yIsDominated) {
//			List<Double> newExternal = new ArrayList<Double>();
//			newExternal.add(f1Val);
//			newExternal.add(f2Val);
//			newExternal.add(weightVector[0]);
//			newExternal.add(weightVector[1]);
//			EP.add(newExternal);
//		}
//	}
//
//	/**
//	 * @param a
//	 * @return true if a dominates b
//	 */
//	private boolean dominates(double[] a, double[] b) {
//		return (a[0] <= b[0] && a[1] <= b[1]) && (a[0] < b[0] || a[1] < b[1]);
//	}

}
