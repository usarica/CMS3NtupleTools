
#ifndef ESCluster_h
#define ESCluster_h

class ESCluster {

	public:

		ESCluster() {}

		ESCluster(float eta, float phi, float energyInPlane1, float energyInPlane2, int usedPlanes)
			: eta_(eta), phi_(phi), 
			energyInPlane1_(energyInPlane1), energyInPlane2_(energyInPlane2), 
			usedPlanes_(usedPlanes) {}

		~ESCluster() {};

		float energy(int plane = 3) {
			if (plane == 3) return energyInPlane1_ + energyInPlane2_;
			if (plane == 1) return energyInPlane1_;
			if (plane == 2) return energyInPlane2_;
			return -1.0;
		}
		int usedPlanes() { return usedPlanes_; }
		float eta() { return eta_; }
		float phi() { return phi_; }

	private:

		float eta_;
		float phi_;
		float energyInPlane1_;
		float energyInPlane2_;
		int usedPlanes_;

};

#endif

