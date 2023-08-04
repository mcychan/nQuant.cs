using System;
using System.Collections.Generic;

namespace nQuant.Master.Ga
{
	public interface Chromosome<T> where T : Chromosome<T>
	{
		public T MakeNewFromPrototype();

		public float Fitness { get; }

		public T Crossover(T mother, int numberOfCrossoverPoints, float crossoverProbability);

		public void Mutation(int mutationSize, float mutationProbability);

		public double[] Objectives { get; }

		public double[] ConvertedObjectives { get; }

		public void ResizeConvertedObjectives(int numObj);

		public Random Random { get; }

		public bool Dominates(T other);

	}
}
