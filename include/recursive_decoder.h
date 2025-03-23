#include "decoder.h"

class recursive_decoder : public soft_decoder {
public:
    recursive_decoder(matrix const&);

private:
    class node {
    protected:
        node(unsigned x, unsigned y, unsigned k1, unsigned k2, unsigned k3, matrix const& Gp);

    public:
        void clear();
        void print(std::ostream& out) const;
        std::string getName() const;

		virtual void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) = 0;
		virtual void downward_pass(std::vector<double>& L) const = 0;

	public:
		const unsigned x, y;
		const unsigned k1, k2, k3;
		std::vector<double> A, B;
		matrix Gp;

	protected:
		inline static const double I = 0;
		virtual void __clear() = 0;
		virtual void __print(std::ostream&) const = 0;
    };

    struct leaf : public node {
		leaf(unsigned x, unsigned y, unsigned k1, matrix const& Gp);

	protected:
		void __print(std::ostream&) const override;
		double F(binvector const& c, std::vector<double> const& p0, std::vector<double> const& p1) const;
        double LLR(double M0, double M1) const;
    };

    struct leaf_no_simplify : public leaf {
		leaf_no_simplify(unsigned x, unsigned y, unsigned k1, matrix const& Gp);

    protected:
		std::vector<std::vector<double>> A1;

		void __clear() override;
        void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) override;
		void downward_pass(std::vector<double>& L) const override;
    };

#ifdef ENABLE_OPT_1
	struct leaf_simplify_1 : public leaf {
		leaf_simplify_1(unsigned x, unsigned y, matrix const& Gp);

    protected:
		double ext_l;

		void __clear() override;
        void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) override;
   		void downward_pass(std::vector<double>& L) const override;
	};
#endif // ENABLE_OPT_1

#ifdef ENABLE_OPT_2
	struct leaf_simplify_2 : public leaf {
		leaf_simplify_2(unsigned x, unsigned y, matrix const& Gp);

    protected:
		void __clear() override;
		void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) override;
		void downward_pass(std::vector<double>& L) const override;
	};
#endif // ENABLE_OPT_2

#ifdef ENABLE_OPT_3
	struct leaf_simplify_3 : public leaf {
		leaf_simplify_3(unsigned x, unsigned y, unsigned k1, matrix const& Gp);

    protected:
		double Phi00, Phi10, Phi01, Phi11;

		void __clear() override;
		void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) override;
		void downward_pass(std::vector<double>& L) const override;
	};
#endif // ENABLE_OPT_3

#ifdef ENABLE_OPT_5
	struct leaf_simplify_5 : public leaf {
		leaf_simplify_5(unsigned x, unsigned y, matrix const& Gp);

    protected:
		void __clear() override;
		void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) override;
		void downward_pass(std::vector<double>& L_out) const override;
	};
#endif // ENABLE_OPT_5

	struct inner : public node {
		inner(unsigned x, unsigned y, unsigned z, node* left, node* right, unsigned k1, unsigned k2, unsigned k3, matrix const& G_hat, matrix const& G_tilda, matrix const& Gp);

    protected:
		unsigned z;
		node *left, *right;
		matrix G_hat, G_tilda;

		void __clear();
		void __print(std::ostream& out) const override;
		void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) override;
		void downward_pass(std::vector<double>& L) const override;
	};

public:
	std::vector<double> decode_soft(std::vector<double> const&) override;
    void printTree(std::string const& fileName) const;

private:
	std::vector<double> p0, p1;
	std::vector<std::vector<unsigned>> _split;
	std::vector<std::vector<unsigned>> _ctor;
	std::vector<binvector> G;
	node* root;

    node* sectionalization();
    node* rec_build_section_tree(unsigned x, unsigned y);
};
