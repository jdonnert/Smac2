enum knot_types {
	KNOT_EMPTY,
	KNOT_START,		// 8bit x, 16bit y, 16bit dxdy, dxddy=0 
	KNOT_MINMAX,		// 3x8bit x, 16bit y0; y1=y2=0 
	KNOT_FULL,		// 3x8bit x, 3x16bit y
	KNOT_STOP,		// 8bit x, 16bit y, 16bit dxdy, dxddy=0
};

struct Knot {			// uncompressed knot
	enum knot_types type;
	int8_t idx;		// knot index
	float P[2];		// knot coordinates
	float Mleft[2];		// left control point 
	float Mright[2];	// right control point 
	float dR;		// holds right derivative of START/STOP
	float dL;		// holds left derivative for START/STOP
};

#pragma pack(push)		// push current packing val to stack
#pragma pack(1)			// tell compiler not to use padding here
struct half_point {		// 40 bit compressed, type HALF_KNOT
	uint8_t x;
	uint16_t y;
	union control_points {
		uint8_t xLR[2];	// extrema have x values only 
		uint8_t xyLR[2];	// extrema have x values only 
	} M;
};

struct full_point {		// 72 bit compressed, type FULL_KNOT
	uint8_t x;
	uint16_t y;
	uint8_t xleft;
	uint8_t xright;
	uint16_t yleft;
	uint16_t yright;
};
#pragma pack(pop)		// pop/restore former packing value

static size_t find_files();
static int uncompress_knots_binary(const char *, struct Knot *);
static float uncompressFloat_16bit(uint16_t), uncompressFloat_8bit(uint8_t);

void Uncompress(const float, const char *, double *);
