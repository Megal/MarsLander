//  Created by Svyatoshenko "Megal" Misha on 2016-04-28

import Foundation

getpid()
//system("lsof -p \(getpid()) >&2")
//system("cat \"\(__FILE__)\" >&2")

public struct StderrOutputStream: OutputStreamType {
	public mutating func write(string: String) { fputs(string, stderr) }
}
public var errStream = StderrOutputStream()

func log(message: String) {	debugPrint("I: "+message, toStream: &errStream) }
@noreturn func fatal(message: String = "E: Some fatal error") { log(message); abort() }

if let inputFile = NSBundle.mainBundle().pathForResource("input", ofType: "txt") {
	freopen(inputFile, "r", stdin)
}

typealias Int2d = (x: Int, y: Int)
func +(a: Int2d, b: Int2d) -> Int2d { return (a.x+b.x, a.y+b.y) }
func -(a: Int2d, b: Int2d) -> Int2d { return (a.x-b.x, a.y-b.y) }

func dist<T where T:Comparable, T:Strideable>(a: T, _ b: T) -> T.Stride {
	if a < b {
		return b - a
	}
	else {
		return a - b
	}
}

func sign<T where T:IntegerLiteralConvertible, T:Comparable>(number: T) -> T {
	if number >= 0 {
		return 1
	} else {
		return -1
	}
}

func clamp(number: Int, range: Range<Int>) -> Int? {
	guard let min = range.minElement(), max = range.maxElement() else {
		return nil
	}
	if number < min {
		return min
	} else if number >= max {
		return max
	} else {
		return number
	}
}

func makeClosedIntervalWithEpsilonMargin(from a: Double, to b: Double, epsilon: Double = 1e-9) -> ClosedInterval<Double>?
{
	if a < b {
		return (a-epsilon ... b+epsilon)
	} else if b < a {
		return (b-epsilon ... a+epsilon)
	} else {
		return nil
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
		let input = input.componentsSeparatedByString(" ").flatMap { Int($0) }
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

extension Range where Element : Comparable {

	init(interval: ClosedInterval<Element>) {
		startIndex = interval.start
		endIndex = interval.end.successor()
	}

	init(interval: HalfOpenInterval<Element>) {
		startIndex = interval.start
		endIndex = interval.end
	}
}

extension MarsLander {

	struct Action {
		let rotate, power: Int
	}

	func clampedAction(action: Action) -> Action {
		let powerClamp = Range(interval: (0...4 as ClosedInterval).clamp(power-1 ... power+1))
		let rotateClamp = Range(interval: (-90...90 as ClosedInterval).clamp(rotate-15 ... rotate+15))

		let clamped = Action(
			rotate: clamp(action.rotate, range: rotateClamp)!,
			power: clamp(action.power, range: powerClamp)!
		)
		return clamped
	}
}

extension MarsLander: CustomStringConvertible {

	static let OneDigitAfterDecimalPoint: (Double) -> String = {
		let formatter = NSNumberFormatter()
		formatter.minimumFractionDigits = 1
		formatter.maximumFractionDigits = 1

		return { (number: Double) -> String in formatter.stringFromNumber(NSNumber(double: number))! }
	}()

	var description: String {
		let d1 = MarsLander.OneDigitAfterDecimalPoint
		return "at:(\(d1(X)), \(d1(Y))) v:(\(d1(hSpeed)), \(d1(vSpeed))) rot:\(rotate) pow:\(power) fuel:\(fuel)"
	}
}

infix operator ~== { associativity left precedence 130 }
//! Equals in acceptable precision
func ~==(left: MarsLander, right: MarsLander) -> Bool {
	guard left.fuel == right.fuel && left.rotate == right.rotate && left.power == right.power else {
		return false
	}

	let epsilon = 1.0 + 1e-9
	guard fabs(left.X - right.X) < epsilon else {
		return false
	}
	guard fabs(left.Y - right.Y) < epsilon else {
		return false
	}
	guard fabs(left.vSpeed - right.vSpeed) < epsilon else {
		return false
	}
	guard fabs(left.hSpeed - right.hSpeed) < epsilon else {
		return false
	}

	return true
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
			let inputs = (readLine()!).componentsSeparatedByString(" ")
			let landX = Int(inputs[0])! // X coordinate of a surface point. (0 to 6999)
			let landY = Int(inputs[1])! // Y coordinate of a surface point. By linking all the points together in a sequential fashion, you form the surface of Mars.
			surface.append((x: landX, y:landY))
		}

		return World(surface: surface)
	}

	init(surface initSurface: [Int2d]) {
		surface = initSurface

		for point in surface {
			maxY = [maxY, point.y].sort(>)[0]
		}

		for (i, point1) in surface.enumerate() where i<surface.count-1 {
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

	func simulate(marsLander lander: MarsLander, action actionBeforeClamp: MarsLander.Action) -> MarsLander {
		let action = lander.clampedAction(actionBeforeClamp)

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
			if let x12 = makeClosedIntervalWithEpsilonMargin(from: x1, to: x2) where x12.contains(lander.X) {
				let t = (lander.X - x1) / (x2 - x1)
				let yt = y1 + t*(y2 - y1)

				return lander.Y > yt
			}
		}

		return true
	}

	func testLanded(marsLander lander: MarsLander) -> Bool {
		let x12 = makeClosedIntervalWithEpsilonMargin(from: target.x-500, to: target.x+500)!
		let y12 = makeClosedIntervalWithEpsilonMargin(from: target.y-40, to: target.y)!
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
			genes.insert((Action(rotate: 0, power: 0), duration: Chromosome.maxTTL), atIndex: 0)
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
		for (currentIndex, sample) in current.enumerate() {
			guard sample.ttl < 0 else { continue }

			evalTTL(&current[currentIndex])
		}
	}

	private mutating func evalTTL(inout chromosome: Chromosome) {
		var evolvingLander = lander
		var ttl = 0
		for (geneIndex, (action: action, duration: duration)) in chromosome.genes.enumerate() {
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
		return Chromosome.Action(rotate: random() % 181 - 90, power: random() % 5)
	}

	mutating func aaa() {
		current.removeAll()
		current.reserveCapacity(populationLimit * 2)
		for _ in 0..<populationLimit {
			let newChromosome = Chromosome(
				genes: [(action: generateRandomAction(), duration: Chromosome.maxTTL)],
				ttl: -1)
			current.append(newChromosome)
		}

		evalTTL()
	}

	mutating func reducePopulation() {
		guard current.count > populationLimit else { return }

		current.sortInPlace{ (a, b) in a.ttl > b.ttl }
		current.removeRange(populationLimit..<current.count)
	}
}

struct RandomDoubleGenerator {
	private init() { srand48(Int(arc4random())) }
	static let singleton = RandomDoubleGenerator()
}

//! Get double in desired interval
extension RandomDoubleGenerator {

	struct Arc4Ranges {
		static let СlosedDenominator: Double = Double(UInt32.max)
		static let HalfOpenDenominator: Double = Double(Int64(UInt32.max) + 1)
	}

	subscript(interval: ClosedInterval<Double>) -> Double {
		let normalized = Double(arc4random()) / RandomDoubleGenerator.Arc4Ranges.СlosedDenominator
		let width = interval.end - interval.start
		let scaled = normalized * width

		return interval.start + scaled
	}

	subscript(interval: HalfOpenInterval<Double>) -> Double {
		let normalized = Double(arc4random()) / RandomDoubleGenerator.Arc4Ranges.HalfOpenDenominator
		let width = interval.end - interval.start
		let scaled = normalized * width

		return interval.start + scaled
	}
}

extension IntervalType {

	public var hashValue: Int {
		if let start = self.start as? NSObject, let end = self.end as? NSObject {
			let halfshift = sizeof(Int)*4
			return start.hashValue ^ ((end.hashValue << halfshift) | (end.hashValue >> halfshift))
		} else {
			return 0
		}
	}
}

extension HalfOpenInterval: Hashable {}
extension ClosedInterval: Hashable {}

//("a"..."b" as ClosedInterval).hashValue
//(0...1 as ClosedInterval).hashValue
//(1e-2...1e+3 as ClosedInterval).hashValue
//
//(0.0..<1.0).hashValue
//(1.0..<2.0).hashValue
//(2.0..<3.0).hashValue
//(0.0..<1.0)==(1.0..<2.0)
//
//var dictionary = [
//	0.0..<1.0 : "Okay",
//	1.0..<2.0 : "Better",
//	2.0..<3.0 : "Perfect"]
//var dict2: Dictionary<HalfOpenInterval<Double>, String> = [
//	0.0..<1.0 : "Okay",
//	1.0..<2.0 : "Better",
//	2.0..<3.0 : "Perfect"]
//
//var dict3: Dictionary<HalfOpenInterval<Double>, String> = [:]
//dict3[0.0..<1.0] = "Meh"
//
//for (range, value) in dict3 {
//	print("\(value) is assign to \(range)")
//}


struct WeightedRandom<T> {

	private var probabilityMap: [HalfOpenInterval<Double> : T] = [:]
	private var weightSum = 0.0
}

enum WeightedRandomError : ErrorType {
	case InvalidArgument
}

extension WeightedRandom {

	mutating func add(value: T, weight: Double) throws {
		guard weight.isNormal && weight > 0 else { throw WeightedRandomError.InvalidArgument }

		let newWeightSum = weightSum + weight
		probabilityMap[weightSum ..< newWeightSum] = value
		weightSum = newWeightSum
	}

	func getRandomObject() -> T? {
		guard !probabilityMap.isEmpty else { return nil }

		let randomizer = RandomDoubleGenerator.singleton
		randomizer[0 ..< weightSum]
		// TODO:
		return nil
	}
}

let world = World.parseFromInput()
var action: MarsLander.Action!
var landerExpected: MarsLander!
for turn in 0..<3000 {
	defer {
		print("\(action.rotate) \(action.power)")
		log("Expect Mars Lander: \(landerExpected)")
	}

	guard feof(stdin) == 0 else {
		break
	}
	let line = readLine()!
	var lander = MarsLander(parseFromInput: line)!
	log("Read   Mars Lander: \(lander)")
	if let landerExpected = landerExpected {
		if lander ~== landerExpected {
			lander = landerExpected
		} else {
			log("Not close to expected, using data from input")
		}
	}

	let testAngle = { (angle: Int) in dist( lander.rotate, angle ) < 90 ? true : false }

	if turn == 0 {
		log("Calculating traectory \((lander.X, lander.Y)) -> \(world.target): ...")
	}

	let altitude = lander.Y - Double(world.target.y)
	// alt*a + (40*40)/2 = v*v/2
	let minVSpeed = 7-Int(sqrt(Double(altitude)*world.a_land*2.0 + 1500))
	log("altitude = \(altitude) minVSpeed = \(minVSpeed)")

	if altitude + lander.vSpeed < 4 - lander.vSpeed {
		log("prepare for hard landing!!!")
		action = MarsLander.Action(rotate: 0, power: 0)
		continue
	}

	let freefall = Int(world.freefallTime(vSpeed: Double(lander.vSpeed), altitude: Double(altitude))*1.2+5)

	// final stage
	switch (lander.hSpeed, lander.vSpeed, altitude) {
	case _ where lander.vSpeed >= 0 : action = MarsLander.Action(rotate: 0, power: 0)
	case _ where lander.vSpeed < -38 : action = MarsLander.Action(rotate: 0, power: 4)
	case _ where lander.hSpeed < -2: action = MarsLander.Action(rotate: -8, power: 3)
	case _ where lander.hSpeed > 2: action = MarsLander.Action(rotate: 8, power: 3)
	case _ where (-30 ..< 0) ~= lander.vSpeed: action = MarsLander.Action(rotate: 0, power: 0)
	case _ where (-38 ..< -30) ~= lander.vSpeed: action = MarsLander.Action(rotate: 0, power: 3)
	default: action = MarsLander.Action(rotate: 0, power: 4)
	}

	landerExpected = world.simulate(marsLander: lander, action: action)
}
