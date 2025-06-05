#ifndef RECURSIVE_DECODER_H
#define RECURSIVE_DECODER_H

#include "decoder.h"
#include "gray_code.h"

namespace short_domain {
	struct recursive_decoder : soft_decoder {
		recursive_decoder(unsigned, matrix const&);

	private:
		class node {
		protected:
			node(unsigned x, unsigned y, unsigned k1, unsigned k2, unsigned k3, matrix const& Gp);

		public:
			void print(std::ostream& out) const;
			std::string getName() const;

			virtual void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) = 0;
			virtual void downward_pass(std::vector<double>& L) = 0;

		public:
			const unsigned x, y;
			const unsigned k1, k2, k3;
			std::vector<double> A, B;
			matrix Gp;
			gray_code traverse_order;

		protected:
			inline static const double I = 0;
			virtual void __print(std::ostream&) const = 0;
		};

		struct leaf : public node {
			leaf(unsigned x, unsigned y, unsigned k1, matrix const& Gp);

		protected:
			void __print(std::ostream&) const override;
			double F(std::size_t const& c, std::vector<double> const& p0, std::vector<double> const& p1) const;
			double LLR(double M0, double M1) const;
		};

		struct leaf_full_code : public leaf {
			leaf_full_code(unsigned x, unsigned y);

			void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) override;
			void downward_pass(std::vector<double>& L) override;
		
		private:
			std::vector<double> L_ans;
		};

		struct leaf_no_simplify : public leaf {
			leaf_no_simplify(unsigned x, unsigned y, unsigned k1, matrix const& Gp);

		protected:
			std::vector<double> P0, P1;
			std::vector<std::vector<double>> A0, A1;

			void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) override;
			void downward_pass(std::vector<double>& L) override;
		};

	#ifdef ENABLE_OPT_1
		struct leaf_simplify_1 : public leaf {
			leaf_simplify_1(unsigned x, unsigned y, matrix const& Gp);

		protected:
			double ext_l;

			void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) override;
			void downward_pass(std::vector<double>&) override;
		};
	#endif // ENABLE_OPT_1

	#ifdef ENABLE_OPT_2
		struct leaf_simplify_2 : public leaf {
			leaf_simplify_2(unsigned x, unsigned y, matrix const& Gp);

		protected:
			void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) override;
			void downward_pass(std::vector<double>&) override;

		private:
			double ext_l;
		};
	#endif // ENABLE_OPT_2

	#ifdef ENABLE_OPT_3
		struct leaf_simplify_3 : public leaf {
			leaf_simplify_3(unsigned x, unsigned y, unsigned k1, matrix const& Gp);

		protected:
			double Phi00, Phi10, Phi01, Phi11;

			void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) override;
			void downward_pass(std::vector<double>&) override;
		};
	#endif // ENABLE_OPT_3

	#ifdef ENABLE_OPT_5
		struct leaf_simplify_5 : public leaf {
			leaf_simplify_5(unsigned x, unsigned y, matrix const& Gp);

		protected:
			void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) override;
			void downward_pass(std::vector<double>&) override;

		private:
			std::vector<double> P0, P1;
		};
	#endif // ENABLE_OPT_5

		struct inner : public node {
			inner(unsigned x, unsigned y, unsigned z, node* left, node* right, unsigned k1, unsigned k2, unsigned k3, matrix const& G_hat, matrix const& G_tilda, matrix const& Gp);

		protected:
			unsigned z;
			node *left, *right;
			matrix G_hat, G_tilda;

			void __print(std::ostream& out) const override;
			void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) override;
			void downward_pass(std::vector<double>&) override;
		};

	public:
		binvector encode(binvector const&) override;
		matrix generate_matrix() const override;
		std::vector<double> decode_soft(std::vector<double> const&) override;
		void printTree(std::string const& fileName) const;

	private:
		matrix G, G_enc;
		std::vector<double> p0, p1;
		std::vector<std::vector<unsigned>> _split;
		std::vector<std::vector<unsigned>> _ctor;
		node* root;

		std::vector<double> L_out;

		void sectionalization();
		node* rec_build_section_tree(unsigned x, unsigned y);
	};
}

#endif // RECURSIVE_DECODER_H