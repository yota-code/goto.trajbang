
Object.defineProperty(Array.prototype, "last", {
	value: function() {
		return this[this.length - 1];
	}
})

function m_copysign(magnitude, sign) {
	return magnitude * ((sign >= 0.0) ? (1.0) : (-1.0)) ;
}