# IsoPeak

A single metabolite will result in multiple peaks in a mass spec run due to isotopologues and adducts yet we normally try to quantify a metabolite associated with a single peak.  By simultaneously searching for multiple peaks of a compound and its adducts, this approach attempted to increase both how accurately we can assign a peak to a specific metabolite and improve the quantification of metabolite abundance by combining the information from multiple peaks.

I abandoned this project because after demonstrating proof of principle, I wasn't that interested in refining the method.  As is, the  approach I came up with is too computationally-intensive, and does not scale-well to large datasets. If anyone is interested in improving this method, feel free.
