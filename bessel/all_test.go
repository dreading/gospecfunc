// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package bessel_test

import (
	. "github.com/dreading/gospecfunc/bessel"
	"math"
	"math/cmplx"
	"testing"
)

func TestBesselI(t *testing.T) {
	testCases := []struct {
		α    float64
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0, 5, 27.23987182360444689454423207588441928247906183222098381517},
		{0, 0.4 + 0.1i, 1.037751775187953907718589364295309545291505675072069744 + 0.02037712677480810300172139986765507486164198939561617240i},
		{0.25, 0.4 + 0.1i, 0.7636677889558093934973812584746228980011733000823678786 + 0.05894567689760461852351886060104048793278243522220770649i},
		{1, 0.4 + 0.1i, 0.2032605049177812708944870481129207354193884289297667737 + 0.05296680165856376304875068492944952721943227874651157176i},
		{0, 333 - 876i, -2.42466995791572885758160970181684018139584605311695e142 - 4.86259747781671556106159806889776193425074468393636e142i},
		{0.5, -1 - 1i, 0.64183847533798587174998976767935220304805492627921709682 - 0.72698064596355457190337225790772836580894088300839807068i},
		{2.5, -1 - 1i, 0.1090512513931498998447281861112676727138680045257758394 + 0.06470001733681579022174741587303933404535864764225205616i},
		{-2.5, -1 - 1i, 0.62942959427995171635310765661865843857908838721710437513 - 0.68414399252960233834658890151304218599337951433338558263i},
		{0, 0, 1},
		{0.5, 0, 0},
	}

	for _, tc := range testCases {
		ζ := I(tc.α, tc.x)
		if veryclose(real(ζ), real(tc.y)) == false {
			t.Fatalf("real(I(%v, %v)): expected %v, got %v", tc.α, tc.x, real(tc.y), real(ζ))
		}
		if veryclose(imag(ζ), imag(tc.y)) == false {
			t.Fatalf("imag(I(%v,%v)): expected %v, got %v", tc.α, tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestBesselIx(t *testing.T) {
	testCases := []struct {
		α    float64
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0, 1, 0.465759607593640436501901529563209998732865911828799168039},
		{0, 0.4 + 0.1i, 0.6956258177175556696817972130134154947747646814992337268 + 0.01365919655776342579472144935485867621944712595959603325i},
		{2.5, 0.4 + 0.1i, 0.00320607871296305168167241107138783625697297503671553856 + 0.00227989345013970917024339680669783780981684359879477301i},
		//{1, 0.4 + 0.1i, 0.2032605049177812708944870481129207354193884289297667737 + 0.05296680165856376304875068492944952721943227874651157176i},
		//{0, 333 - 876i, -2.42466995791572885758160970181684018139584605311695e142 - 4.86259747781671556106159806889776193425074468393636e142i},
		//{ 0.5, -1-1i, 0.64183847533798587174998976767935220304805492627921709682 - 0.72698064596355457190337225790772836580894088300839807068i},
		//{2.5, -1-1i, 0.1090512513931498998447281861112676727138680045257758394 + 0.06470001733681579022174741587303933404535864764225205616i},
		{-2.5, -1 - 1i, 0.23155420740047630584325533438599797859706097185198654487 - 0.25168250965258951856814168520313162092554516265352057022i},
		{0, 0, 1},
		{0.5, 0, 0},
	}

	for _, tc := range testCases {
		ζ := Ix(tc.α, tc.x)
		if close(real(ζ), real(tc.y)) == false {
			t.Fatalf("real(Ix(%v, %v)): expected %v, got %v", tc.α, tc.x, real(tc.y), real(ζ))
		}
		if close(imag(ζ), imag(tc.y)) == false {
			t.Fatalf("imag(Ix(%v,%v)): expected %v, got %v", tc.α, tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestBesselJ(t *testing.T) {
	testCases := []struct {
		α    float64
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0, 1, 0.765197686557966551449717526102663220909274289755325241861},
		{0, 0.4 + 0.1i, 0.9627513455152844523412642515338186980389471552835888233 - 0.01962711629303271170774431370757696151493347606833454781i},
	}

	for _, tc := range testCases {
		ζ := J(tc.α, tc.x)
		if veryclose(real(ζ), real(tc.y)) == false {
			t.Fatalf("real(J(%v, %v)): expected %v, got %v", tc.α, tc.x, real(tc.y), real(ζ))
		}
		if veryclose(imag(ζ), imag(tc.y)) == false {
			t.Fatalf("imag(J(%v,%v)): expected %v, got %v", tc.α, tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestBesselJx(t *testing.T) {
	testCases := []struct {
		α    float64
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0, 1, 0.765197686557966551449717526102663220909274289755325241861},
		{0, 0.4 + 0.1i, 0.8711334416866959908400496970749705623358862660709028054 - 0.01775934923007923297551626999816498584343265028733313589i},
	}

	for _, tc := range testCases {
		ζ := Jx(tc.α, tc.x)
		if veryclose(real(ζ), real(tc.y)) == false {
			t.Fatalf("real(J(%v, %v)): expected %v, got %v", tc.α, tc.x, real(tc.y), real(ζ))
		}
		if veryclose(imag(ζ), imag(tc.y)) == false {
			t.Fatalf("imag(J(%v,%v)): expected %v, got %v", tc.α, tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestBesselK(t *testing.T) {
	testCases := []struct {
		α    float64
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0, 1, 0.421024438240708333335627379212609036136219748226660472298},
		{0, 0.4 + 0.1i, 1.0826035097235081563303065645729426351352961828491217647 - 0.21324459634740555852107854608689805665579172058695006929i},
		{2.5, -1 - 1i, 0.81740838955020351274340488272166961021804282493424837140 - 1.1762814200405309189093806024582991851908631533909975900i},
		{-2.5, -1 - 1i, 0.81740838955020351274340488272166961021804282493424837140 - 1.1762814200405309189093806024582991851908631533909975900i},
		{2.5, -1 + 1i, 0.81740838955020351274340488272166961021804282493424837140 + 1.1762814200405309189093806024582991851908631533909975900i},
		{2.5, 1 + 1i, -0.9730203208880581731806799234864631370588047705275199210 - 1.160002999791696943340858523335033177428742515466007164i},
		{25, 1 + 1i, 1.243920310213260987851496108498392250170999248421594e27 - 1.296860880863851282936708596501597343689035562826055e27i},
	}

	for _, tc := range testCases {
		ζ := K(tc.α, tc.x)
		if close(real(ζ), real(tc.y)) == false {
			t.Fatalf("real(K(%v, %v)): expected %v, got %v", tc.α, tc.x, real(tc.y), real(ζ))
		}
		if close(imag(ζ), imag(tc.y)) == false {
			t.Fatalf("imag(K(%v,%v)): expected %v, got %v", tc.α, tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestBesselY(t *testing.T) {
	testCases := []struct {
		α    float64
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0, 1, 0.088256964215676957982926766023515162827817523090675546711},
		{0, 0.4 + 0.1i, -0.5873828735984329430094809243367886365339042810690287965 + 0.1750446663657107385528975496319911274395118490119615635i},
	}

	for _, tc := range testCases {
		ζ := Y(tc.α, tc.x)
		if close(real(ζ), real(tc.y)) == false {
			t.Fatalf("real(K(%v, %v)): expected %v, got %v", tc.α, tc.x, real(tc.y), real(ζ))
		}
		if close(imag(ζ), imag(tc.y)) == false {
			t.Fatalf("imag(K(%v,%v)): expected %v, got %v", tc.α, tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestBesselYx(t *testing.T) {
	testCases := []struct {
		α    float64
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0, 1, 0.088256964215676957982926766023515162827817523090675546711},
		{0, 0.4 + 0.1i, -0.5314860027453484704174295757203952414847036250133517715 + 0.1583869639553156798942662170390978238327603923087624910i},
	}

	for _, tc := range testCases {
		ζ := Yx(tc.α, tc.x)
		if close(real(ζ), real(tc.y)) == false {
			t.Fatalf("real(K(%v, %v)): expected %v, got %v", tc.α, tc.x, real(tc.y), real(ζ))
		}
		if close(imag(ζ), imag(tc.y)) == false {
			t.Fatalf("imag(K(%v,%v)): expected %v, got %v", tc.α, tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestBesselH1(t *testing.T) {
	testCases := []struct {
		α    float64
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0, 1, 0.7651976865579665514497175261026632209092742897553252418 + 0.08825696421567695798292676602351516282781752309067554671i},
		{-1, 1, -0.4400505857449335159596822037189149131273723019927652511 + 0.7812128213002887165471500000479648205499063907164446078i},
		{-0.75, 1, 0.04470111581450463105545772617920751765494009770171649803 + 0.8347550483584058646827795439479346446696466075007092333i},
		{-1, 0.5 + 0.75i, 0.31538103527285508005364790245087868485498240454260042539 + 0.41218087155579729290725721358335683019274777412959566968i},
		{2, 0.75 - 0.75i, 1.0798261992096498876295467929919458008519826441628797726 - 0.54090933484017059182686757461203629371839519480775444538i},
		{2, 0.75 + 0.75i, -1.053476390787392933122805327012746855287254952784286092 - 0.2605861026050603823088046945658155986159626229320885120i},
		{1.25, 75 + 75i, -1.9247251946295874088544114732839933719832641099026e-34 + 7.9990631314254957714076851577823797996920648281887e-35i},
		{0.5, -1 - 1i, 0.323099244779188375242885442046800758630027408670702 - 1.794951464930692261137043373938578831872923335613907i},
		{2.5, 0.75 - 0.75i, 1.9460186316088613793436289685117901696545662575971048949 + 0.31385846950840913303716441049258209193668876274758270688i},
		{0, 1e-19, 0.99999999999999999999999999999999999999750000000000000000 - 27.9253570525269413784402661820401607251961941715595049764i},
	}

	for _, tc := range testCases {
		ζ := H1(tc.α, tc.x)
		if close(real(ζ), real(tc.y)) == false {
			t.Fatalf("real(H1(%v, %v)): expected %v, got %v", tc.α, tc.x, real(tc.y), real(ζ))
		}
		if close(imag(ζ), imag(tc.y)) == false {
			t.Fatalf("imag(H1(%v,%v)): expected %v, got %v", tc.α, tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestBesselH1x(t *testing.T) {
	testCases := []struct {
		α    float64
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0, 1, 0.48770374908695631836308738133058745712258693077093166555 - 0.59620620960600407144981219216668717558529366315691036762i},
		{2, 1, -1.326918900901948524727452213317969542394493645112173208 - 0.9885555673427596851350488773161207595226003330752383517i},
		{2, 0.75 - 0.75i, 0.19905150622945094148807297673110791707463556268369924020 - 0.53463803588129573648193110763794857626578041708802159850i},
		{2, 0.75 - 0.75i, 0.19905150622945094148807297673110791707463556268369924020 - 0.53463803588129573648193110763794857626578041708802159850i},
		{-2, 1, -1.326918900901948524727452213317969542394493645112173208 - 0.9885555673427596851350488773161207595226003330752383517i},
	}

	for _, tc := range testCases {
		ζ := H1x(tc.α, tc.x)
		if close(real(ζ), real(tc.y)) == false {
			t.Fatalf("real(H1(%v, %v)): expected %v, got %v", tc.α, tc.x, real(tc.y), real(ζ))
		}
		if close(imag(ζ), imag(tc.y)) == false {
			t.Fatalf("imag(H1(%v,%v)): expected %v, got %v", tc.α, tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestBesselH1Inf(t *testing.T) {
	testCases := []struct {
		α float64
		x complex128
	}{
		{0, 0},
	}

	for _, tc := range testCases {
		ζ := H1(tc.α, tc.x)
		if cmplx.IsInf(ζ) == false {
			t.Fatalf("H1(%v,%v): expected +Inf, got %v", tc.α, tc.x, ζ)
		}

	}
}

func TestBesselH2(t *testing.T) {
	testCases := []struct {
		α    float64
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0, 1, 0.7651976865579665514497175261026632209092742897553252418 - 0.08825696421567695798292676602351516282781752309067554671i},
		{0.75, 1, 0.55865249320489174775136094526622318790147389214947110819 + 0.62186941744297463828809818067447936925140957188674420518i},
		{-0.75, 1, 0.04470111581450463105545772617920751765494009770171649803 - 0.8347550483584058646827795439479346446696466075007092333i},
	}

	for _, tc := range testCases {
		ζ := H2(tc.α, tc.x)
		if close(real(ζ), real(tc.y)) == false {
			t.Fatalf("real(H(%v, %v)): expected %v, got %v", tc.α, tc.x, real(tc.y), real(ζ))
		}
		if close(imag(ζ), imag(tc.y)) == false {
			t.Fatalf("imag(H(%v,%v)): expected %v, got %v", tc.α, tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestBesselH2x(t *testing.T) {
	testCases := []struct {
		α    float64
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0, 1, 0.48770374908695631836308738133058745712258693077093166555 + 0.59620620960600407144981219216668717558529366315691036762i},
		{0.75, 1, -0.2214438408600644965517736888448038467127370196949914935 + 0.8060873438158229137295900975781045807643143928147235726i},
		{-0.75, 1, 0.72657426856496667374843507246248076381806979708592493717 - 0.41340538551667411791769533165611069490472126141079781951i},
		//{-0.75, 1,  0.82512630137524102887529822062070573012895112676781594289 - 0.13409238342919102684393053416118248699111082458006453873i},
	}

	for _, tc := range testCases {
		ζ := H2x(tc.α, tc.x)
		if close(real(ζ), real(tc.y)) == false {
			t.Fatalf("real(H(%v, %v)): expected %v, got %v", tc.α, tc.x, real(tc.y), real(ζ))
		}
		if close(imag(ζ), imag(tc.y)) == false {
			t.Fatalf("imag(H(%v,%v)): expected %v, got %v", tc.α, tc.x, imag(tc.y), imag(ζ))
		}
	}
}
func TestAiryAi(t *testing.T) {
	testCases := []struct {
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{5, 1.08344428136074417349865025033459804795777834796889391e-04},
		{5i, 29.9014823980071663793368297194577996107966955412306646481 + 21.6778315987836365558366270522814009526154445079208768325i},
		{-0.8 - 5i, -40.502439645308446470999918851421417821301812631645554975 - 117.58608694901610862936384033209723226387843162614560540i},
		{-80 - 50i, 1.28264717210434580202856048898366187480222379931511e196 + 4.86799742198478384026489758810414265077624185873133e195i},
		{69 + 27.345345i, -4.1680951604712214279525229780753604917334277179166e-158 - 3.1646991120150008455176064001093550612164673982528e-158i},
	}

	for _, tc := range testCases {
		ζ := Ai(tc.x)
		if soclose(real(ζ), real(tc.y), 1e-13) == false {
			t.Fatalf("real(AiryAi(%v)): expected %v, got %v", tc.x, real(tc.y), real(ζ))
		}
		if soclose(imag(ζ), imag(tc.y), 1e-13) == false {
			t.Fatalf("imag(AiryAi(%v)): expected %v, got %v", tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestAiryAix(t *testing.T) {
	testCases := []struct {
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{3 + 4i, 0.1835755300356563207855609335333982182804779386314715657987711225267727422186395137060866322176677229 - 0.0416176487393802789393930914530936600764046614258631476023947559642313271327799016541303441791855652i},
	}

	for _, tc := range testCases {
		ζ := Aix(tc.x)
		if soclose(real(ζ), real(tc.y), 1e-13) == false {
			t.Fatalf("real(AiryAix(%v)): expected %v, got %v", tc.x, real(tc.y), real(ζ))
		}
		if soclose(imag(ζ), imag(tc.y), 1e-13) == false {
			t.Fatalf("imag(AiryAix(%v)): expected %v, got %v", tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestAiryAid(t *testing.T) {
	testCases := []struct {
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{3 + 4i, -0.0752099611959030290360245211620050520074612674725256824 + 0.0823640771555377950900025995378362045029329536404160740i},
	}

	for _, tc := range testCases {
		ζ := Aid(tc.x)
		if soclose(real(ζ), real(tc.y), 1e-13) == false {
			t.Fatalf("real(AiryAid(%v)): expected %v, got %v", tc.x, real(tc.y), real(ζ))
		}
		if soclose(imag(ζ), imag(tc.y), 1e-13) == false {
			t.Fatalf("imag(AiryAid(%v)): expected %v, got %v", tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestAiryAidx(t *testing.T) {
	testCases := []struct {
		x, y complex128
	}{
		// extended precision values computed using Mathematica (using the principal branch of the logarithm for complex exponentiation)
		{3 + 4i, -0.412990905801422564636582340314435196054059892109907071 - 0.0920837077431263403171214472631110596105237477987462376i},
	}

	for _, tc := range testCases {
		ζ := Aidx(tc.x)
		if soclose(real(ζ), real(tc.y), 1e-13) == false {
			t.Fatalf("real(AiryAid(%v)): expected %v, got %v", tc.x, real(tc.y), real(ζ))
		}
		if soclose(imag(ζ), imag(tc.y), 1e-13) == false {
			t.Fatalf("imag(AiryAid(%v)): expected %v, got %v", tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestAiryBi(t *testing.T) {
	testCases := []struct {
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{5, 657.7920441711711824410805788744438785563124063292868570873},
	}

	for _, tc := range testCases {
		ζ := Bi(tc.x)
		if soclose(real(ζ), real(tc.y), 1e-13) == false {
			t.Fatalf("real(AiryBi(%v)): expected %v, got %v", tc.x, real(tc.y), real(ζ))
		}
		if soclose(imag(ζ), imag(tc.y), 1e-13) == false {
			t.Fatalf("imag(AiryBi(%v)): expected %v, got %v", tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestAiryBix(t *testing.T) {
	testCases := []struct {
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{3 + 4i, 0.27319149262040082252117234505365461539534889567258960373 + 0.27713977915811109029797070575343074622129648765459620679i},
	}

	for _, tc := range testCases {
		ζ := Bix(tc.x)
		if soclose(real(ζ), real(tc.y), 1e-13) == false {
			t.Fatalf("real(AiryBix(%v)): expected %v, got %v", tc.x, real(tc.y), real(ζ))
		}
		if soclose(imag(ζ), imag(tc.y), 1e-13) == false {
			t.Fatalf("imag(AiryBix(%v)): expected %v, got %v", tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestAiryBid(t *testing.T) {
	testCases := []struct {
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{3 + 4i, 0.78788923789635748275993347147037942948472756149505740013 + 2.9998668872583759608026286983418251783375077031064431087i},
	}

	for _, tc := range testCases {
		ζ := Bid(tc.x)
		if soclose(real(ζ), real(tc.y), 1e-13) == false {
			t.Fatalf("real(AiryBid(%v)): expected %v, got %v", tc.x, real(tc.y), real(ζ))
		}
		if soclose(imag(ζ), imag(tc.y), 1e-13) == false {
			t.Fatalf("imag(AiryBid(%v)): expected %v, got %v", tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestAiryBidx(t *testing.T) {
	testCases := []struct {
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{3 + 4i, 0.20768534826166084976217984558838878267020909768365666904 + 0.79075632620944147532503042962473032501930744949927623196i},
	}

	for _, tc := range testCases {
		ζ := Bidx(tc.x)
		if soclose(real(ζ), real(tc.y), 1e-13) == false {
			t.Fatalf("real(AiryBidx(%v)): expected %v, got %v", tc.x, real(tc.y), real(ζ))
		}
		if soclose(imag(ζ), imag(tc.y), 1e-13) == false {
			t.Fatalf("imag(AiryBidx(%v)): expected %v, got %v", tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestGi(t *testing.T) {
	testCases := []struct {
		num, den, res float64
	}{
		{-1.0e0, 512.0e0, 0.20468308070040542435e0},
		{-1.0e0, 8.0e0, 0.18374662832557904078e0},
		{-1.0e0, 1.0e0, -0.11667221729601528265e0},
		{-4.0e0, 1.0e0, 0.31466934902729557596e0},
		{-8.0e0, 1.0e0, -0.37089040722426257729e0},
		{-33.0e0, 4.0e0, -0.25293059772424019694e0},
		{-9.0e0, 1.0e0, 0.28967410658692701936e0},
		{-10.0e0, 1.0e0, -0.34644836492634090590e0},
		{-11.0e0, 1.0e0, 0.28076035913873049496e0},
		{-13.0e0, 1.0e0, 0.21814994508094865815e0},
		{1.0e0, 512.0e0, 0.20526679000810503329e0},
		{1.0e0, 8.0e0, 0.22123695363784773258e0},
		{1.0e0, 1.0e0, 0.23521843981043793760e0},
		{4.0e0, 1.0e0, 0.82834303363768729338e-1},
		{7.0e0, 1.0e0, 0.45757385490989281893e-1},
		{29.0e0, 4.0e0, 0.44150012014605159922e-1},
		{8.0e0, 1.0e0, 0.39951133719508907541e-1},
		{9.0e0, 1.0e0, 0.35467706833949671483e-1},
		{10.0e0, 1.0e0, 0.31896005100679587981e-1},
		{12.0e0, 1.0e0, 0.26556892713512410405e-1},
		{-1.0e40, 1.0e0, 0},
		{1.0e-16, 1.0e0, 0.20497554248200024505},
		{-1.0e-16, 1.0e0, 0.20497554248200024505},
		{1.0e20, 1.0e0, 0},
		{1.0e10, 1.0e0, 3.183098861837907e-11},
		{1.0e308, 1.0e0, 0},
		{-1.0e8, 1.0e0, -0.0009916070219422016},
		{0e0, 1e0, 0.204975542482000245050307456364537851198242729549532168346},
	}

	for _, tc := range testCases {
		ζ := Gi(tc.num / tc.den)
		if close(ζ, tc.res) == false {
			t.Fatalf("Gi(%v): expected %v, got %v", tc.num/tc.den, tc.res, ζ)
		}

	}
}

func TestHi(t *testing.T) {
	testCases := []struct {
		num, den, res float64
	}{
		{-1.0e0, 512.0e0, 0.40936798278458884024e0},
		{-1.0e0, 8.0e0, 0.37495291608048868619e0},
		{-1.0e0, 1.0e0, 0.22066960679295989454e0},
		{-4.0e0, 1.0e0, 0.77565356679703713590e-1},
		{-8.0e0, 1.0e0, 0.39638826473124717315e-1},
		{-33.0e0, 4.0e0, 0.38450072575004151871e-1},
		{-9.0e0, 1.0e0, 0.35273216868317898556e-1},
		{-10.0e0, 1.0e0, 0.31768535282502272742e-1},
		{-11.0e0, 1.0e0, 0.28894408288051391369e-1},
		{-13.0e0, 1.0e0, 0.24463284011678541180e-1},
		{1.0e0, 512.0e0, 0.41053540139998941517e0},
		{1.0e0, 8.0e0, 0.44993502381204990817e0},
		{1.0e0, 1.0e0, 0.97220515514243332184e0},
		{4.0e0, 1.0e0, 0.83764237105104371193e2},
		{7.0e0, 1.0e0, 0.80327744952044756016e5},
		{29.0e0, 4.0e0, 0.15514138847749108298e6},
		{8.0e0, 1.0e0, 0.11995859641733262114e7},
		{9.0e0, 1.0e0, 0.21472868855967642259e8},
		{10.0e0, 1.0e0, 0.45564115351632913590e9},
		{12.0e0, 1.0e0, 0.32980722582904761929e12},
		{-1.0e40, 1.0e0, 0},
		{1.0e-16, 1.0e0, 0.40995108496400049010},
		{-1.0e-16, 1.0e0, 0.40995108496400049010},
		{-1.0e8, 1.0e0, 3.183098861837907e-09},
		{-1.0e308, 1.0e0, 0},
		{0e0, 1e0, 0.409951084964000490100614912729075702396485459099064336693},
		{100e0, 1e0, 6.0412239966702013990052650701823592789077760491387386e288},
	}

	for _, tc := range testCases {
		ζ := Hi(tc.num / tc.den)
		if close(ζ, tc.res) == false {
			t.Fatalf("Hi(%v): expected %v, got %v", tc.num/tc.den, tc.res, ζ)
		}

	}
}

func TestHiInf(t *testing.T) {
	testCases := []struct {
		val float64
	}{
		{104e0},
		{1.0e308},
	}

	for _, tc := range testCases {
		ζ := Hi(tc.val)
		if math.IsInf(ζ, 1) == false {
			t.Fatalf("Hi(%v): expected +Inf, got %v", tc.val, ζ)
		}

	}
}

// The floating point comparison tests are copied from from math/all_test.go.
func tolerance(a, b, e float64) bool {
	// Multiplying by e here can underflow denormal values to zero.
	// Check a==b so that at least if a and b are small and identical
	// we say they match.
	if a == b {
		return true
	}
	d := a - b
	if d < 0 {
		d = -d
	}
	// note: b is correct (expected) value, a is actual value.
	// make error tolerance a fraction of b, not a.
	if b != 0 {
		e = e * b
		if e < 0 {
			e = -e
		}
	}
	return d < e
}

func close(a, b float64) bool      { return tolerance(a, b, 1e-14) }
func veryclose(a, b float64) bool  { return tolerance(a, b, 5e-16) }
func soclose(a, b, e float64) bool { return tolerance(a, b, e) }
