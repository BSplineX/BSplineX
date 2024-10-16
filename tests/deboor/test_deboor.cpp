// Third-party includes
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// BSplineX includes
#include "deboor/deboor.hpp"

using namespace Catch::Matchers;
using namespace bsplinex;
using namespace bsplinex::deboor;

TEST_CASE(
    "deboor::Deboor<T, C, BoundaryCondition::OPEN, EXT> deboor{knots::Knots<T, "
    "C, BoundaryCondition::OPEN, EXT> knots, "
    "control_points::ControlPoints<T, BoundaryCondition::OPEN> control_points, "
    "degree}",
    "[deboor]"
)
{
  size_t degree{3};

  std::vector<double> t_data_vec{0.1, 1.3, 2.2, 2.2, 4.9, 6.3, 6.3, 6.3, 13.2};
  knots::Data<double, Curve::NON_UNIFORM> t_data{t_data_vec};
  knots::Knots<
      double,
      Curve::NON_UNIFORM,
      BoundaryCondition::OPEN,
      Extrapolation::NONE>
      knots{t_data, degree};

  std::vector<double> c_data_vec{0.1, 1.3, 2.2, 4.9, 13.2};
  control_points::Data<double> c_data{c_data_vec};
  control_points::ControlPoints<double, BoundaryCondition::OPEN> control_points{
      c_data, degree
  };

  DeBoor<
      double,
      Curve::NON_UNIFORM,
      BoundaryCondition::OPEN,
      Extrapolation::NONE>
      deboor{knots, control_points, degree};

  // To generate this reference data use the `reference.py` script in this
  // folder
  // clang-format off
  std::vector<double> x_values{2.2, 2.3000000000000003, 2.4000000000000004, 2.5000000000000004, 2.6000000000000005, 2.7000000000000006, 2.8000000000000007, 2.900000000000001, 3.000000000000001, 3.100000000000001, 3.200000000000001, 3.300000000000001, 3.4000000000000012, 3.5000000000000013, 3.6000000000000014, 3.7000000000000015, 3.8000000000000016, 3.9000000000000017, 4.000000000000002, 4.100000000000001, 4.200000000000002, 4.3000000000000025, 4.400000000000002, 4.500000000000002, 4.600000000000002, 4.700000000000003, 4.8000000000000025, 4.900000000000002, 5.000000000000003, 5.100000000000003, 5.200000000000003, 5.3000000000000025, 5.400000000000003, 5.5000000000000036, 5.600000000000003, 5.700000000000003, 5.800000000000003, 5.900000000000004, 6.0000000000000036, 6.100000000000003, 6.200000000000004};
  std::vector<double> y_values{0.4, 0.4987905929445725, 0.5953834608104188, 0.6901102371457323, 0.7833025554987059, 0.8752920494175335, 0.9664103524504085, 1.0569890981455239, 1.1473599200510731, 1.2378544517152497, 1.3288043266862466, 1.420541178512258, 1.5133966407414765, 1.6077023469220952, 1.7037899306023085, 1.8019910253303089, 1.9026372646542904, 2.0060602821224456, 2.112591711282968, 2.2225631856840518, 2.3363063388738903, 2.4541528044006755, 2.5764342158126015, 2.7034822066578617, 2.83562841048465, 2.9732044608411603, 3.1165419912755823, 3.2659726353361123, 3.4243850625148546, 3.604896086079547, 3.8231795552418366, 4.094909319213373, 4.435759227205808, 4.861403128430789, 5.387514872099959, 6.029768307424969, 6.8038372836174785, 7.725395649889126, 8.810117255451555, 10.073675949516419, 11.53174558129538};
  // clang-format on

  SECTION("deboor.deboor(...)")
  {
    for (size_t i{0}; i < x_values.size(); i++)
    {
      REQUIRE_THAT(
          deboor.deboor(knots.find(x_values.at(i)), x_values.at(i)),
          WithinRel(y_values.at(i))
      );
    }
  }
}

TEST_CASE(
    "deboor::Deboor<T, C, BoundaryCondition::CLAMPED, EXT> "
    "deboor{knots::Knots<T, C, BoundaryCondition::CLAMPED, EXT> knots, "
    "control_points::ControlPoints<T, BoundaryCondition::CLAMPED> "
    "control_points, degree}",
    "[deboor]"
)
{
  size_t degree{3};

  std::vector<double> t_data_vec{0.1, 1.3, 2.2, 2.2, 4.9, 6.3, 6.3, 6.3, 13.2};
  knots::Data<double, Curve::NON_UNIFORM> t_data{t_data_vec};
  knots::Knots<
      double,
      Curve::NON_UNIFORM,
      BoundaryCondition::CLAMPED,
      Extrapolation::NONE>
      knots{t_data, degree};

  std::vector<double> c_data_vec{
      0.1, 1.3, 2.2, 2.5, 3.2, 4.3, 4.9, 5.6, 5.4, 0.3, 13.2
  };
  control_points::Data<double> c_data{c_data_vec};
  control_points::ControlPoints<double, BoundaryCondition::CLAMPED>
      control_points{c_data, degree};

  DeBoor<
      double,
      Curve::NON_UNIFORM,
      BoundaryCondition::CLAMPED,
      Extrapolation::NONE>
      deboor{knots, control_points, degree};

  // To generate this reference data use the `reference.py` script in this
  // folder
  // clang-format off
  std::vector<double> x_values{0.1, 0.2, 0.30000000000000004, 0.4, 0.5, 0.6000000000000001, 0.7000000000000001, 0.8, 0.9, 1.0, 1.1, 1.2000000000000002, 1.3, 1.4000000000000001, 1.5, 1.6, 1.7000000000000002, 1.8, 1.9000000000000001, 2.0, 2.1, 2.2, 2.3000000000000003, 2.4000000000000004, 2.5, 2.6, 2.7, 2.8000000000000003, 2.9000000000000004, 3.0, 3.1, 3.2, 3.3000000000000003, 3.4000000000000004, 3.5, 3.6, 3.7, 3.8000000000000003, 3.9000000000000004, 4.0, 4.1000000000000005, 4.2, 4.3, 4.4, 4.5, 4.6000000000000005, 4.7, 4.800000000000001, 4.9, 5.0, 5.1000000000000005, 5.2, 5.300000000000001, 5.4, 5.5, 5.6000000000000005, 5.7, 5.800000000000001, 5.9, 6.0, 6.1000000000000005, 6.2, 6.300000000000001, 6.4, 6.5, 6.6000000000000005, 6.7, 6.800000000000001, 6.9, 7.0, 7.1000000000000005, 7.2, 7.300000000000001, 7.4, 7.5, 7.6000000000000005, 7.7, 7.800000000000001, 7.9, 8.0, 8.1, 8.200000000000001, 8.3, 8.4, 8.5, 8.6, 8.700000000000001, 8.8, 8.9, 9.0, 9.1, 9.200000000000001, 9.3, 9.4, 9.5, 9.600000000000001, 9.700000000000001, 9.8, 9.9, 10.0, 10.100000000000001, 10.200000000000001, 10.3, 10.4, 10.5, 10.600000000000001, 10.700000000000001, 10.8, 10.9, 11.0, 11.100000000000001, 11.200000000000001, 11.3, 11.4, 11.5, 11.600000000000001, 11.700000000000001, 11.8, 11.9, 12.0, 12.100000000000001, 12.200000000000001, 12.3, 12.4, 12.5, 12.600000000000001, 12.700000000000001, 12.8, 12.9, 13.0, 13.100000000000001};
  std::vector<double> y_values{0.1, 0.38599773242630386, 0.6451247165532882, 0.8790816326530615, 1.0895691609977327, 1.2782879818594108, 1.4469387755102046, 1.5972222222222225, 1.7308390022675746, 1.849489795918368, 1.954875283446712, 2.0486961451247168, 2.13265306122449, 2.2083781249125165, 2.2772290809327846, 2.3404950869236583, 2.399465300523502, 2.455428879370679, 2.509674981103552, 2.5634927633604874, 2.6181713837798495, 2.675, 2.7341327228722783, 2.7947796390456796, 2.8568138698550487, 2.920108536635234, 2.984536760721082, 3.049971663447438, 3.116286366149151, 3.1833539901610672, 3.2510476568180313, 3.319240487454893, 3.387805603406498, 3.4566161260076926, 3.525545176593323, 3.594465876498238, 3.663251347057283, 3.731774709605305, 3.7999090854771507, 3.8675275960076667, 3.934503362531702, 4.0007095063841005, 4.066019148899709, 4.130305411413377, 4.193441415259949, 4.255300281774272, 4.315755132291193, 4.3746790881455615, 4.4319452706722196, 4.487649631536136, 4.542779443722759, 4.598544810547657, 4.6561558353263965, 4.716822621374545, 4.781755272007673, 4.852163890541345, 4.929258580291129, 5.014249444572595, 5.108346586701308, 5.212760109992836, 5.328700117762751, 5.457376713326616, 5.599999999999999, 5.588286470081489, 5.5708160202612405, 5.54800690392044, 5.520277374440273, 5.488045685201927, 5.451730089586586, 5.411748840975436, 5.368520192749664, 5.322462398290456, 5.273993710978998, 5.223532384196474, 5.171496671324073, 5.118304825742976, 5.064375100834375, 5.010125749979451, 4.955975026559393, 4.902341183955385, 4.849642475548615, 4.798297154720266, 4.748723474851525, 4.70133968932358, 4.656564051517614, 4.614814814814814, 4.576510232596366, 4.542068558243458, 4.511908045137272, 4.486446946658996, 4.466103516189815, 4.4512960071109156, 4.442442672803485, 4.439961766648707, 4.444271542027768, 4.455790252321855, 4.474936150912152, 4.502127491179847, 4.537782526506123, 4.582319510272169, 4.63615669585917, 4.699712336648312, 4.773404686020779, 4.857651997357759, 4.952872524040438, 5.0594845194500015, 5.177906236967634, 5.308555929974522, 5.451851851851853, 5.608212255980811, 5.7780553957425855, 5.961799524518356, 6.159862895689313, 6.37266376263664, 6.600620378741526, 6.844150997385161, 7.103673871948718, 7.37960725581339, 7.6723694023603635, 7.982378564970824, 8.310052997025965, 8.655810951906956, 9.02007068299499, 9.403250443671258, 9.80576848731694, 10.228043067313234, 10.670492437041306, 11.133534849882354, 11.617588559217563, 12.123071818428114, 12.650402880895207};
  // clang-format on

  SECTION("deboor.deboor(...)")
  {
    for (size_t i{0}; i < x_values.size(); i++)
    {
      REQUIRE_THAT(
          deboor.deboor(knots.find(x_values.at(i)), x_values.at(i)),
          WithinRel(y_values.at(i))
      );
    }
  }
}
