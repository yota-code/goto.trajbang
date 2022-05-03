

class TrajBang3 {

	constructor() {
		
	}

	compute(jm, am, a0, s0, sg) {

		var k = (0 <= (sg - s0)) ? (1.0) : (-1.0);
		if ( (a0**2) / (2*k*jm) + s0 == sg ) {
			return [
				[-k],
				[(a0) / (k*jm)],
				'A'
			];
		}

		for ( let k of [-1, 1] ) {
			var m = (a0**2 / 2) + k*jm*(sg - s0);
			if ( 0 <= m ) {
				var q = Math.sqrt(m);
				var d1 = (-a0 + k*q) / (k*jm);
				var d2 = (k*q) / (k*jm);
				if ( 0 <= d1 && 0 <= d2 && q <= am ) {
					return [
						[k, -k],
						[ d1, d2 ],
						'B' + ((0 < k) ? ('+') : ('-'))
					];
				}
			}
		}
		
		for ( let k of [-1, 1] ) {
			var c0 = (0 <= k*am - a0) ? (1.0) : (-1.0);
			var d0 = (k*am - a0) / (c0*jm);
			var s01 = d0*(a0+k*am) / 2;

			var c2 = (0 <= -k*am) ? (1.0) : (-1.0);
			var d2 = (-k*am)/(c2*jm);
			var s23 = d2*(k*am) / 2;

			var s12 = sg - s01 - s23 - s0;
			var d1 = s12 / (k*am);

			if ( 0.0 <= d1 ) {
				return [
					[c0, 0, c2],
					[d0, d1, d2],
					'C' + ((0 < k) ? ('+') : ('-'))
				];
			}
		}

		return [ [0, 0, 0], [1, 1, 1], 'Z' ];

	}

	update() {

		var jm = parseFloat(document.forms.tb3_form.elements.tb3_jm_in.value);
		var am = parseFloat(document.forms.tb3_form.elements.tb3_am_in.value);
		var a0 = parseFloat(document.forms.tb3_form.elements.tb3_a0_in.value);
		var s0 = parseFloat(document.forms.tb3_form.elements.tb3_s0_in.value);
		var sg = parseFloat(document.forms.tb3_form.elements.tb3_sg_in.value);
	
		document.forms.tb3_form.elements.tb3_jm_out.value = `${jm}`;
		document.forms.tb3_form.elements.tb3_am_out.value = `${am}`;
		document.forms.tb3_form.elements.tb3_a0_out.value = `${a0}`;
		document.forms.tb3_form.elements.tb3_s0_out.value = `${s0}`;
		document.forms.tb3_form.elements.tb3_sg_out.value = `${sg}`;
	
		var [cmd, dur, bch] = this.compute(jm, am, a0, s0, sg);
	
		console.log(dur);
		console.log(cmd)
	
		document.forms.tb3_form.elements.tb3_bch_out.value = `${bch}`;
	
		var result_lst = new Array();
		var [Tp, Ap, Sp, Pp] = [0.0, a0, s0, 0.0];
		for ( let i=0 ; i<cmd.length ; i++ ) {
			var c = cmd[i];
			var d = dur[i];
			var result = this.unroll(d, Tp, c*jm, Ap, Sp, Pp, 32);
			var [Tp, Ap, Sp, Pp] = [
				result[0].last(),
				result[2].last(),
				result[3].last(),
				result[4].last()
			];
			result_lst.push(result);
		}
	
		var j_plot = new Array();
		var a_plot = new Array();
		var s_plot = new Array();
		var p_plot = new Array();
	
		for ( let i=0 ; i<result_lst.length ; i++ ) {
	
			var j =	{ x: result_lst[i][0], y: result_lst[i][1], mode: 'lines' };
			j_plot.push(j);
			var a = { x: result_lst[i][0], y: result_lst[i][2], mode: 'lines' };
			a_plot.push(a);
			var s = { x: result_lst[i][0], y: result_lst[i][3], mode: 'lines' };
			s_plot.push(s);
			var p = { x: result_lst[i][0], y: result_lst[i][4], mode: 'lines' };
			p_plot.push(p);
	
		}
	
		var Time = result_lst.last()[0].last();
		document.forms.tb3_form.elements.tb3_Time_out.value = `${Time.toFixed(3)}`;
	
		var Pg = result_lst.last()[4].last();
		document.forms.tb3_form.elements.tb3_Pg_out.value = `${Pg.toFixed(3)}`;
	
		var layout = {'showlegend': false};
	
		Plotly.newPlot('j_plot', j_plot, layout);
		Plotly.newPlot('a_plot', a_plot, layout);
		Plotly.newPlot('s_plot', s_plot, layout);
		Plotly.newPlot('p_plot', p_plot, layout);
	
	}

	unroll(d, Tp, Jc, Ap, Sp, Pp, n) {

		var t_lst = new Array();
	
		var j_lst = new Array();
		var a_lst = new Array();
		var s_lst = new Array();
		var p_lst = new Array();
	
		for ( let i=0 ; i<n+1 ; i++ ) {
			var t = (d * i)/n;
			t_lst.push( t + Tp );
			j_lst.push( Jc );
			a_lst.push( Ap + Jc*t );
			s_lst.push( Sp + Ap*t + (1.0/2.0)*Jc*t**2 );
			p_lst.push( Pp + Sp*t + (1.0/2.0)*Ap*t**2 + (1.0/6.0)*Jc*t**3 );
		}
	
		return [t_lst, j_lst, a_lst, s_lst, p_lst];
	
	}

}


