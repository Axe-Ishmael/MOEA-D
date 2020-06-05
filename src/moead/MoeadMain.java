package moead;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

public class MoeadMain {

	private static final int ITERATIONS = 100;
	private static final int POPULATION_SIZE = 30;
	private static final int NEIGHBOR_SIZE = 10;
	private static final int GENOME_SIZE = 10;

	public static void main(String[] args) {
		if (POPULATION_SIZE <= NEIGHBOR_SIZE) {
			System.out.println("Population size (" + POPULATION_SIZE
					+ ") must be greater/equal to neighbor size ("
					+ NEIGHBOR_SIZE + "), aborting...");
			System.exit(1);
		}

		Moead moead = new Moead();
		List<List<Double>> solution = moead.moead(ITERATIONS, POPULATION_SIZE, NEIGHBOR_SIZE, GENOME_SIZE);
		Printer.printSolution(solution);
//		saveToFile("/desktop/moead-sol.txt", solution);
		saveToFile("C:\\Users\\Administrator\\Desktop", solution);
	}

	private static void saveToFile(String filename, List<List<Double>> solution) {
		try {
			BufferedWriter br = new BufferedWriter(new FileWriter(filename));
			for (List<Double> l : solution) {
				br.write(l.get(0) + " " + l.get(1) + "\n");
			}
			br.flush();
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
