using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using HaploGrepSharp;
using HaploGrepSharp.NewSearchMethods;
using System.Linq;

namespace HaploGrepSharp.TreeUtilities
{
	/// <summary>
	/// A wrapper class that manages the addition and removal of polymorphisms as we move down the phylotree, it 
	/// is somewhat complicated because of SNPs and insertions/deletions appearing at the same position, and so requiring multiple
	/// matching hits
	/// </summary>
	public class PolymorphismCollection : IEnumerable<Polymorphism>
	{
		private Dictionary<int, MutationalHistory> polymorphisms = new Dictionary<int, MutationalHistory> ();

		public void AddRange (PolymorphismCollection otherPolys)
		{
			foreach (var v in otherPolys)
				Add (v);
		}

		public void FilterSites (PolymorphismFilter filterToApply)
		{
			var goodPolys = filterToApply.FilterPolys (this).Select (x => x.Position).ToDictionary (z => z, u => u);
			List<int> toRemove = this.Select (x => x.position).Where (y => !goodPolys.ContainsKey (y)).ToList ();
			foreach (var v in toRemove) {
				polymorphisms.Remove (v);		
			}		
		}

		public void Add (Polymorphism poly)
		{
			if (!polymorphisms.ContainsKey (poly.position)) {
				polymorphisms [poly.position] = new MutationalHistory (poly);
			} else {
				polymorphisms [poly.position].AddPoly (poly);
			}
		}

		public IEnumerator<Polymorphism> GetEnumerator ()
		{
			return polymorphisms.Values.Where (z => z.CurrentState != null).SelectMany (x => x.CurrentState).GetEnumerator ();
		}

		System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator ()
		{
			return GetEnumerator ();
		}

		/// <summary>
		/// A class that holds the mutational history of all polymorphisms that occurred at this site from the rCRS on the way to this value
		/// </summary>
		public class MutationalHistory
		{
			public static List<List<Polymorphism>> UnresolvedSeriesOfChanges = new List<List<Polymorphism>> ();
			/// <summary>
			/// The current state relative to the rCRS, if a mutation is undone by a later mutation, this is an empty list
			/// 
			/// </summary>
			public List<Polymorphism> CurrentState;
			/// <summary>
			/// The ordered list of all changes that appeared at this position
			/// </summary>
			public List<Polymorphism> History;

			public MutationalHistory (Polymorphism p)
			{
				History = new List<Polymorphism> () { p };
				CurrentState = new List<Polymorphism> () { p };
			}

			public void AddPoly (Polymorphism p)
			{
				History.Add (p);
				//This is a very backward way of representing this, as far as I can tell from the source code
				//if 152C on the tree is followed by 152C!, then this means the 152C mutation should simply be cleared
				if (p.isBackMutation) {
					//make sure we have one to remove
					var toDrop = CurrentState.Where (x => x.position == p.position && x.mutation == p.mutation).ToList ();
					if (toDrop.Count != 1) {
						throw new HaploGrepException ("Cannot back mutate when no mutation appeared!");
					}
					CurrentState.Remove (toDrop [0]);
				} else {                    
					var c1Count = CurrentState.Count == 1;
					var c1SNP = MutationAssigner.MutationIsBasePair (CurrentState [0].mutation);
					//simple SNP replacement
					if (c1Count && c1SNP) {
						CurrentState.Clear ();
						CurrentState.Add (p);
					}
                    //insertions are appended on, and once any deletion happens we are shit out of luck
                    else if (p.mutation == Mutations.INS && c1SNP && c1Count) {
						CurrentState.Add (p);
					} else if (p.mutation == Mutations.INS) {
						//Some spot have repeated inserts of the same base, e.g. 
						if (CurrentState.Any (z => z.mutation != Mutations.INS) || CurrentState.Any (z => z.InsertedPolys.Length > 1)) {
							throw new HaploGrepException ("Can't add insertions over multiple complex backgrounds");
						} else {
							//TODO: 455.1T and 455.2T are simply indicative of a two bp insertion (at the 1 position and the 2 position, this is totally confusiing).
							//Verify that all insertions are the same, this happens in some cases such as 455.1T being followed by a 455.2T, note I know realize this 
							var set = new HashSet<string> (CurrentState.Select (x => x.insertedPolys));
							set.Add (p.insertedPolys);
							if (set.Count > 1) {
								throw new Exception ("Cannot add a new type of insertion");
							}
							CurrentState.Add (p); 
						}
					}
                    //Special exception for 1719.1G followed by 1719A, which is another odd case
                    else if (p.position == 1719 && p.ToString () == "1719A") {
						if (CurrentState.Count != 1 || CurrentState [0].ToString () != "1719.1G") {
							throw new HaploGrepException ("1719 posiiton exception");
						}
						CurrentState.Add (p);
						CurrentState.Reverse ();
					}
                    //this can happen if the base of interest is not the reference or last mutation.                
                    else {
						if (!UnresolvedSeriesOfChanges.Contains (History)) {
							UnresolvedSeriesOfChanges.Add (History);
						}
						throw new HaploGrepException ("Tried to change a polymorphism at a location where a mutation already occurred without it being asimple thing!");
					}
				}
			}
		}
	}
}
