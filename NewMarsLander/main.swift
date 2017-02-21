//  Created by Svyatoshenko "Megal" Misha on 2016-04-28

import Foundation

// MARK: - Logging stuff
public struct StderrOutputStream: TextOutputStream {
	public mutating func write(_ string: String) { fputs(string, stderr) }
}
public var errStream = StderrOutputStream()
func log(_ message: String) {	print(message, to: &errStream) }
func fatal(_ message: String = "Fatal error!") -> Never  { log(message); abort() }

// Local Tests
if let inputFile = Bundle.main.path(forResource: "input", ofType: "txt") {
	freopen(inputFile, "r", stdin)
}

// Cartesian 2d
typealias Int2d = (x: Int, y: Int)
func +(a: Int2d, b: Int2d) -> Int2d { return (a.x+b.x, a.y+b.y) }
func -(a: Int2d, b: Int2d) -> Int2d { return (a.x-b.x, a.y-b.y) }


// MARK: - Range extensions

extension Range where Bound: FloatingPoint {

	/// Adds an ε-neighboorhood to create a safe zone for calculation errors
	func extended(by epsilon: Bound = Bound.ulpOfOne.squareRoot()) -> Range {

		return (self.lowerBound.nextDown-epsilon ..< self.upperBound.nextUp+epsilon)
	}
}

extension ClosedRange where Bound: FloatingPoint {

	/// Adds an ε-neighboorhood to create a safe zone for calculation errors
	func extended(by epsilon: Bound = Bound.ulpOfOne.squareRoot()) -> ClosedRange {

		return (self.lowerBound.nextDown-epsilon ... self.upperBound.nextUp+epsilon)
	}
}

extension CountableClosedRange {

	/// Clamps a value in range
	func clamp(_ value: Bound) -> Bound {

		if value <= lowerBound {
			return lowerBound
		}
		else if value >= upperBound {
			return upperBound
		}
		else {
			return value
		}
	}
}

struct MarsLander {

	var X: Double
	var Y: Double
	var hSpeed: Double
	var vSpeed: Double
	var fuel: Int
	var rotate: Int
	var power: Int
}

extension MarsLander {

	init?(parseFromInput input: String) {
		let input = input.components(separatedBy: " ").flatMap { Int($0) }
		guard input.count == 7 else {
			return nil
		}

		self.init(
			X: Double(input[0]),
			Y: Double(input[1]),
			hSpeed: Double(input[2]),
			vSpeed: Double(input[3]),
			fuel: input[4],
			rotate: input[5],
			power: input[6]
		)
	}
}



extension MarsLander {

	struct Action {
		let rotate, power: Int
	}

	func clamped(action: Action) -> Action {

		let rotateClamp = (-90...90).clamped(to: rotate-15 ... rotate+15)
		let powerClamp = (0...4).clamped(to: power-1 ... power+1)

		return Action(
			rotate: rotateClamp.clamp(action.rotate),
			power: powerClamp.clamp(action.power)
		)
	}
}


func -(_ a: MarsLander, _ b: MarsLander) -> Double {

	let dvx = abs(a.hSpeed - b.hSpeed)
	let dvy = abs(a.vSpeed - b.vSpeed)
	let dx = abs(a.X - b.X)
	let dy = abs(a.Y - b.Y)

	return dx + dy + dvx + dvy
}


extension MarsLander: CustomStringConvertible {

	static let doubleFormatter: (Double) -> String = {
		let formatter = NumberFormatter()
		formatter.minimumFractionDigits = 1
		formatter.maximumFractionDigits = 1

		return { (number: Double) -> String in formatter.string(from: NSNumber(value: number as Double))! }
	}()

	var description: String {
		let d1 = MarsLander.doubleFormatter
		return "at:(\(d1(X)), \(d1(Y))) v:(\(d1(hSpeed)), \(d1(vSpeed))) rot:\(rotate) pow:\(power) fuel:\(fuel)"
	}
}

struct World {
	var surface: [Int2d]
	var maxY: Int = 0
	var target: (x: Double, y: Double) = (0, 0)
	let gravity = 3.711
	let a_land = 4.0 - 3.711
	let a_8 = 0.55669240384026

	static func parseFromInput() -> World {
		let surfaceN = Int(readLine()!)! // the number of points used to draw the surface of Mars.
		var surface = [Int2d]()
		for _ in 0..<surfaceN {
			let inputs = (readLine()!).components(separatedBy: " ")
			let landX = Int(inputs[0])! // X coordinate of a surface point. (0 to 6999)
			let landY = Int(inputs[1])! // Y coordinate of a surface point. By linking all the points together in a sequential fashion, you form the surface of Mars.
			surface.append((x: landX, y:landY))
		}

		return World(surface: surface)
	}

	init(surface initSurface: [Int2d]) {
		surface = initSurface

		for point in surface {
			maxY = [maxY, point.y].sorted(by: >)[0]
		}

		for (i, point1) in surface.enumerated() where i<surface.count-1 {
			let point2 = surface[i+1]
			if point1.y == point2.y {
				target = (x: Double(point1.x + point2.x)/2, y: Double(point1.y))
				break
			}
		}
	}
}

extension World {

	func freefallTime(vSpeed v1: Double, altitude h: Double) -> Double {
		let v2 = sqrt( v1*v1 + 2*gravity*h )
		let t = (v2 - v1) / gravity

		return t
	}

	func simulate(marsLander lander: MarsLander, action beforeClamp: MarsLander.Action) -> MarsLander {

		let action = lander.clamped(action: beforeClamp)

		let trueAngle = Double(90+action.rotate) * M_PI/180.0
		let (dvx, dvy) = (
			Double(action.power) * cos(trueAngle),
			Double(action.power) * sin(trueAngle) - gravity
		)

		var next = lander
		next.fuel -= action.power
		next.hSpeed += dvx
		next.vSpeed += dvy
		next.X += (lander.hSpeed + next.hSpeed) / 2
		next.Y += (lander.vSpeed + next.vSpeed) / 2
		next.power = action.power
		next.rotate = action.rotate

		return next
	}

	func testSafeZone(marsLander lander: MarsLander) -> Bool {
		guard (0.0..<7000.0) ~= lander.X else {
			return false
		}
		guard (0.0..<3000.0) ~= lander.Y else {
			return false
		}

		for i in 0..<surface.count-1 {
			let (x1, x2) = (Double(surface[i].x), Double(surface[i+1].x))
			let (y1, y2) = (Double(surface[i].y), Double(surface[i+1].y))
			let x12 = (x1...x2).extended()
			if x12 ~= lander.X {
				let t = (lander.X - x1) / (x2 - x1)
				let yt = y1 + t*(y2 - y1)

				return lander.Y > yt
			}
		}

		return true
	}

	func testLanded(marsLander lander: MarsLander) -> Bool {

		let x12 = (target.x-500 ... target.x+500).extended()
		let y12 = (target.y-40 ... target.y).extended()

		guard x12 ~= lander.X && y12 ~= lander.Y else {
			return false
		}
		guard lander.vSpeed > -40 else {
			return false
		}
		guard fabs(lander.hSpeed) < 20 else {
			return false
		}
		guard lander.rotate == 0 else {
			return false
		}

		return true
	}
}

// MARK: - Random helpers

func xorshift128plus(seed0 : UInt64, seed1 : UInt64) -> () -> UInt64 {
	var s0 = seed0
	var s1 = seed1
	if s0 == 0 && s1 == 0 {
		s1 =  1 // The state must be seeded so that it is not everywhere zero.
	}

	return {
		var x = s0
		let y = s1
		s0 = y
		x ^= x << 23
		x ^= x >> 17
		x ^= y
		x ^= y >> 26
		s1 = x
		return s0 &+ s1
	}

}

struct Random {

	let generator = xorshift128plus(seed0: 232323232323, seed1: 123987345675)

	func bounded(to max: UInt64) -> UInt64 {
		var u: UInt64 = 0
		let b: UInt64 = (u &- max) % max
		repeat {
			u = generator()
		} while u < b
		return u % max
	}

	/// Random value for `Int` in arbitrary closed range, uniformally distributed
	subscript(range: CountableClosedRange<Int>) -> Int {
		let bound = range.upperBound.toIntMax() - range.lowerBound.toIntMax() + 1
		let x = range.lowerBound + Int(bounded(to: UInt64(bound)))

		guard range.contains(x) else { fatal("out of range") }
		return x
	}

	/// Random value for `Double` in arbitrary closed range
	subscript(range: ClosedRange<Double>) -> Double {
		let step = (range.upperBound - range.lowerBound) / Double(UInt64.max)

		let value = range.lowerBound + step*Double(generator())
		guard range.contains(value) else { fatal("out of range") }

		return value
	}

	/// Random value for `Double` in arbitrary half-open range
	subscript(range: Range<Double>) -> Double {
		let step = (range.upperBound - range.lowerBound) / (1.0 + Double(UInt64.max))

		let value = range.lowerBound + step*Double(generator())
		guard range.contains(value) else { fatal("out of range") }

		return value
	}

}

let random = Random()

struct Chromosome {

	typealias Action = MarsLander.Action
	typealias GeneElement = (action: Action, duration: Int)

	var genes = [Chromosome.neutralGene]
	var ttl = -1

	static let neutralGene: GeneElement = (Action(rotate: 0, power: 0), duration: Chromosome.maxTTL)
	static let maxTTL = 3000
}

extension Chromosome {

	mutating func incrementAge() {
		genes[0].duration -= 1
		if genes[0].duration < 1 {
			genes.removeFirst()
		}
		if genes.isEmpty {
			genes.insert((Action(rotate: 0, power: 0), duration: Chromosome.maxTTL), at: 0)
			ttl = -1
		} else {
			ttl -= 1
		}
	}
}

struct Generation {
	var current: [Chromosome] = []
	let populationLimit = 100
	let crossingoverRate = 2
	let mutationRate = 3
	let world: World
	var lander: MarsLander

	init(world: World, lander: MarsLander) {
		self.world = world
		self.lander = lander
		current.append(Chromosome(genes: [Chromosome.neutralGene], ttl: -1))
	}
}

extension Generation {

	mutating func evalTTL() {
		for (currentIndex, sample) in current.enumerated() {
			guard sample.ttl < 0 else { continue }

			evalTTL(&current[currentIndex])
		}
	}

	fileprivate mutating func evalTTL(_ chromosome: inout Chromosome) {
		var evolvingLander = lander
		var ttl = 0
		for (geneIndex, (action: action, duration: duration)) in chromosome.genes.enumerated() {
			for turnSameAction in 0..<duration {
				ttl += 1
				let nextLander = world.simulate(marsLander: evolvingLander, action: action)
				defer {
					evolvingLander = nextLander
				}

				if world.testSafeZone(marsLander: nextLander) { continue }

				let isLanded = world.testLanded(marsLander: nextLander)
				if isLanded {
					chromosome.ttl = Chromosome.maxTTL
					return
				} else {
					let alternativeLander = world.simulate(marsLander: evolvingLander, action: Chromosome.neutralGene.action)
					if world.testLanded(marsLander: alternativeLander) {
						// TODO: cut tail and add neutral gene
						assert(false)
					} else {
						// TODO: cut tail and maybe? add neutral gene
						chromosome.ttl = ttl
						return
					}
				}
			}
		}
		assert(false)
	}
}

extension Generation {

	func generateRandomAction() -> Chromosome.Action {
		return Chromosome.Action(rotate: random[-90...90], power: random[0...4])
	}

	mutating func populateToLimitWithRandom() {
		current.reserveCapacity(populationLimit * 2)
		for _ in current.count..<populationLimit {
			let newChromosome = Chromosome(
				genes: [(action: generateRandomAction(), duration: Chromosome.maxTTL)],
				ttl: -1)
			current.append(newChromosome)
		}

		evalTTL()
	}

	mutating func reducePopulation() {
		reducePopulation(populationLimit)
	}

	mutating func reducePopulation(_ limit: Int) {
		guard current.count > limit else { return }

		current.sort{ (a, b) in a.ttl > b.ttl }
		current.removeSubrange(limit..<current.count)
	}

	func bestChomosome() -> Chromosome { return current.first! }

	mutating func nextTurn(marsLander next: MarsLander) {
		lander = next
		for i in 0..<current.count {
			current[i].incrementAge()
		}
	}
}

//struct RandomDoubleGenerator {
//	private init() { srand48(Int(arc4random())) }
//	static let singleton = RandomDoubleGenerator()
//}
//
////! Get double in desired interval
//extension RandomDoubleGenerator {
//
//	struct Arc4Ranges {
//		static let СlosedDenominator: Double = Double(UInt32.max)
//		static let HalfOpenDenominator: Double = Double(Int64(UInt32.max) + 1)
//	}
//
//	subscript(interval: ClosedInterval<Double>) -> Double {
//		let normalized = Double(arc4random()) / RandomDoubleGenerator.Arc4Ranges.СlosedDenominator
//		let width = interval.end - interval.start
//		let scaled = normalized * width
//
//		return interval.start + scaled
//	}
//
//	subscript(interval: HalfOpenInterval<Double>) -> Double {
//		let normalized = Double(arc4random()) / RandomDoubleGenerator.Arc4Ranges.HalfOpenDenominator
//		let width = interval.end - interval.start
//		let scaled = normalized * width
//
//		return interval.start + scaled
//	}
//}

//extension IntervalType {
//
//	public var hashValue: Int {
//		if let start = self.start as? NSObject, let end = self.end as? NSObject {
//			let halfshift = MemoryLayout<Int>.size*4
//			return start.hashValue ^ ((end.hashValue << halfshift) | (end.hashValue >> halfshift))
//			//infix operator .--> {
//			//associativity left
//			//precedence 152
//			//}
//			////! Apply operation using dot syntax
//			//public func .--> <U, V>(arg: U, transform: (U) -> V ) -> V {
//			//    return transform(arg)
//			//}
//		} else {
//			return 0
//		}
//	}
//}
//
//extension Range: Hashable {}
//extension ClosedRange: Hashable {}

//
////("a"..."b" as ClosedInterval).hashValue
////(0...1 as ClosedInterval).hashValue
////(1e-2...1e+3 as ClosedInterval).hashValue
////
////(0.0..<1.0).hashValue
////(1.0..<2.0).hashValue
////(2.0..<3.0).hashValue
////(0.0..<1.0)==(1.0..<2.0)
////
////var dictionary = [
////	0.0..<1.0 : "Okay",
////	1.0..<2.0 : "Better",
////	2.0..<3.0 : "Perfect"]
////var dict2: Dictionary<HalfOpenInterval<Double>, String> = [
////	0.0..<1.0 : "Okay",
////	1.0..<2.0 : "Better",
////	2.0..<3.0 : "Perfect"]
////
////var dict3: Dictionary<HalfOpenInterval<Double>, String> = [:]
////dict3[0.0..<1.0] = "Meh"
////
////for (range, value) in dict3 {
////	print("\(value) is assign to \(range)")
////}
//
//
////struct WeightedRandom<T> {
////
////	private var probabilityMap: [HalfOpenInterval<Double> : T] = [:]
////	private var weightSum = 0.0
////}
////
////enum WeightedRandomError : ErrorType {
////	case InvalidArgument
////}
////
////extension WeightedRandom {
////
////	mutating func add(value: T, weight: Double) throws {
////		guard weight.isNormal && weight > 0 else { throw WeightedRandomError.InvalidArgument }
////
////		let newWeightSum = weightSum + weight
////		probabilityMap[weightSum ..< newWeightSum] = value
////		weightSum = newWeightSum
////	}
////
////	func getRandomObject() -> T? {
////		guard !probabilityMap.isEmpty else { return nil }
////
////		let randomizer = RandomDoubleGenerator.singleton
////		randomizer[0 ..< weightSum]
////		// TODO:
////		return nil
////	}
////}

let world = World.parseFromInput()
var action: MarsLander.Action!
var landerExpected: MarsLander!
var generation: Generation!
for turn in 0..<3000 {
	defer {
		print("\(action.rotate) \(action.power)")
		log("Expect Mars Lander: \(landerExpected!)")
	}

	guard feof(stdin) == 0 else {
		break
	}
	let line = readLine()!
	var lander = MarsLander(parseFromInput: line)!
	log("Read   Mars Lander: \(lander)")
	if let landerExpected = landerExpected {
		if lander - landerExpected < 1.0 {
			lander = landerExpected
		} else {
			log("Not close to expected, using data from input")
		}
	}

	let testAngle = { (angle: Int) in abs(lander.rotate - angle) < 90 ? true : false }

	if turn == 0 {
		log("Calculating traectory \((lander.X, lander.Y)) -> \(world.target): ...")
		generation = Generation(world: world, lander: lander)
		for _ in 0..<30 {
			generation.reducePopulation(generation.populationLimit/2)
			generation.populateToLimitWithRandom()
			generation.evalTTL()
		}
	} else {
		generation.nextTurn(marsLander: lander)
		generation.evalTTL()
		for _ in 0..<20 {
			generation.reducePopulation(generation.populationLimit/2)
			generation.populateToLimitWithRandom()
			generation.evalTTL()
		}
	}

	action = generation.bestChomosome().genes[0].action

	landerExpected = world.simulate(marsLander: lander, action: action)
}
