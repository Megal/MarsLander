//  Created by Svyatoshenko "Megal" Misha on 2016-04-28

import Foundation

// MARK: - Dot Transform operation

precedencegroup DotOperationPrecedence {
	higherThan: MultiplicationPrecedence
	assignment: true
}
infix operator .--> : DotOperationPrecedence
public func .--> <U, V>(arg: U, transform: (U) -> V ) -> V {
	return transform(arg)
}

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

// MARK: - Cartesian 2d
typealias Int2d = (x: Int, y: Int)
func +(a: Int2d, b: Int2d) -> Int2d { return (a.x+b.x, a.y+b.y) }
func -(a: Int2d, b: Int2d) -> Int2d { return (a.x-b.x, a.y-b.y) }

typealias Double2d = (x: Double, y: Double)
func +(a: Double2d, b: Double2d) -> Double2d { return (a.x+b.x, a.y+b.y) }
func -(a: Double2d, b: Double2d) -> Double2d { return (a.x-b.x, a.y-b.y) }


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

	struct Action: Hashable {
		let rotate, power: Int

		// MARK: - Hashable
		var hashValue: Int {
			return (90 + rotate) + power << 8
		}

		// MARK: - Equatable
		static func ==(lhs: Action, rhs: Action) -> Bool {
			return lhs.rotate == rhs.rotate && lhs.power == rhs.power
		}
	}

	func clamped(action: Action) -> Action {

		let rotateClamp = (-90...90).clamped(to: rotate-15 ... rotate+15)
		let powerClamp = (power-1 ... power+1).clamped(to: 0...fuel).clamped(to: 0...4)

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

//MARK: - DNA

struct Chromosome {

	typealias Action = MarsLander.Action
	typealias GeneElement = (action: Action, duration: Int)

	var genes = [Chromosome.neutralGene]
	var ttl = -1
	var blackBox: MarsLander?

	static let neutralGene: GeneElement = (Action(rotate: 0, power: 0), Chromosome.maxTTL)
	static let maxTTL = 300
}

extension Chromosome {

	var totalLength: Int {

		return genes.reduce(0) { $0 + $1.duration }
	}

	mutating func incrementAge() {
		genes[0].duration -= 1
		if genes[0].duration < 1 {
			genes.removeFirst()
		}
		if genes.isEmpty {
			genes.insert(Chromosome.neutralGene, at: 0)
			ttl = -1
		} else {
			ttl -= 1
			genes[genes.count-1].duration += 1
		}
	}

	func prefix(n: Int) -> Chromosome {
		guard n > 0 else { return Chromosome(genes: [Chromosome.neutralGene], ttl: -1, blackBox: nil) }

		var copyLength = 0
		var copyGenes: [GeneElement] = []
		for gene in genes {
			if copyLength + gene.duration >= n {
				copyGenes.append(GeneElement(action: gene.action, duration: n-copyLength))
				break;
			} else {
				copyGenes.append(gene)
				copyLength += gene.duration
			}
		}

		var prefix = Chromosome(genes: copyGenes, ttl: -1, blackBox: nil)
		prefix.appendNeutralTail()
		return prefix
	}

	func suffix(from: Int) -> Chromosome {
		var skipped = 0
		var copyGenes: [GeneElement] = []
		for gene in genes {
			if skipped + gene.duration < from {
				skipped += gene.duration
			} else {
				let appendLength = skipped + gene.duration - from
				copyGenes.append(GeneElement(action: gene.action, duration: appendLength))
				skipped = from
			}
		}

		if copyGenes.count == 0 { copyGenes.append(Chromosome.neutralGene) }
		return Chromosome(genes: copyGenes, ttl: -1, blackBox: nil)
	}

	static func combine(head: Chromosome, tail: Chromosome, at cutPoint: Int) -> Chromosome {
		var copied: [GeneElement] = []
		var headCopiedLength = 0
		for gene in head.genes {
			if gene.duration + headCopiedLength < cutPoint {
				copied.append(gene)
				headCopiedLength += gene.duration
			} else {
				let cutGene = (action: gene.action, duration: cutPoint - headCopiedLength)
				copied.append(cutGene)
				break;
			}
		}

		var tailSkippedLength = 0
		for gene in tail.genes {
			if tailSkippedLength >= cutPoint {
				copied.append(gene)
			} else if tailSkippedLength + gene.duration > cutPoint {
				let cutLength = tailSkippedLength + gene.duration - cutPoint
				let copingGene = (action: gene.action, duration: cutLength)
				copied.append(copingGene)
				tailSkippedLength += gene.duration
			} else {
				tailSkippedLength += gene.duration
			}
		}

		let combined = Chromosome(genes: copied, ttl: -1, blackBox: nil)
		if combined.totalLength != maxTTL {
			log("chromosome has bad length = \(combined.totalLength)")
		}
		return combined
	}

	mutating func appendNeutralTail() {
		let total = totalLength
		guard total < Chromosome.maxTTL else { return }

		let neutralTail = (action: Chromosome.neutralGene.action, duration: Chromosome.maxTTL-totalLength)
		if genes.isEmpty {
			genes.append(neutralTail)
		}
		else if genes.last!.action == neutralTail.action {
			genes[genes.count-1].duration += neutralTail.duration
		}
		else {
			genes.append(neutralTail)
		}
	}
}


struct Generation {
	var current: [Chromosome] = {
		var current: [Chromosome] = []
		current.reserveCapacity(50)

		return current
	}()
	var fitnessScore: [Double] = []

	let populationLimit = 20
	let crossingoverRate = 2
	let mutationRate = 3
	let world: World
	var lander: MarsLander
	var fitnessScored = false {
		didSet { sorted = false	}
	}

	var sortedIndexes: [Int] = []
	var sorted = false

	var lastEvaluatedPath: [Double2d] = []
	///0 is acceptable, >0 and more is worse
	typealias ErrorFn = (MarsLander) -> Double
	var fitnessFunc: [(fn: ErrorFn, weight: Double)] = []

	init(world: World, lander: MarsLander) {
		self.world = world
		self.lander = lander
		current.append(Chromosome(genes: [Chromosome.neutralGene], ttl: -1, blackBox: nil))
	}
}

extension Generation {

	mutating func evalTTL() {
		guard fitnessScored == false else { return }

		for (currentIndex, sample) in current.enumerated() {
			guard sample.ttl < 0 else { continue }

			evalTTL(&current[currentIndex])
		}

		evalFitness()
	}

	mutating func evalFitness() {
		guard !fitnessScored else { return }
		defer {
			fitnessScored = true
		}

		guard case let N = current.count, N > 0 else { return }
		fitnessScore = [Double](repeating: 1.0, count: N)

		for (fn, denormWeight) in fitnessFunc {
			let weight = denormWeight / Double(N)
			let errors = current.map { $0.blackBox!.-->fn }
			let indexes = (0..<N).map { (i: Int($0), e: errors[$0]) }
			let worstFirst = indexes.sorted { $0.e > $1.e }

			var worse = 0
			for i in 1..<N {
				let prev = worstFirst[i-1]
				let cur = worstFirst[i]

				if prev.e > 1e-9 + cur.e {
					worse = i
				}

				fitnessScore[cur.i] += weight * Double(worse)
			}
		}
	}

	mutating func sort() {
		guard fitnessScored else { assert(false, "you should call evalFitness() before this method"); return }
		guard !sorted else { return }
		defer {
			sorted = true
		}

		guard case let N = current.count, N > 0 else { return }
		sortedIndexes = (0..<N).map { Int($0) }
			.sorted { (i, j) in
				if fitnessScore[i] == fitnessScore[j] {
					return current[i].ttl > current[j].ttl
				} else {
					return fitnessScore[i] > fitnessScore[j]
				}
			}
	}

	mutating func evalBest() {
		evalTTL()
		sort()
		evalTTL(&current[sortedIndexes[0]])
	}

	mutating func evalSuboptimal(place n: Int ) {
		evalTTL()
		sort()
		evalTTL(&current[sortedIndexes[n]])
	}

	private mutating func evalTTL(_ chromosome: inout Chromosome) {
		let oldTTL = chromosome.ttl
		defer {
			assert(chromosome.blackBox != nil, "blackbox is necessary")
			if chromosome.ttl != oldTTL || oldTTL < 0 {
				fitnessScored = false
			}
		}

		lastEvaluatedPath.removeAll()
		var evolvingLander = lander; lastEvaluatedPath.append((x: evolvingLander.X, y: evolvingLander.Y))
		var ttl = 0
		for (geneIndex, (action: action, duration: duration)) in chromosome.genes.enumerated() {
			for turnSameAction in 0..<duration {
				ttl += 1
				let nextLander = world.simulate(marsLander: evolvingLander, action: action)
				defer {
					lastEvaluatedPath.append((x: evolvingLander.X, y: evolvingLander.Y))
					evolvingLander = nextLander
				}

				if world.testSafeZone(marsLander: nextLander) { continue }

				let isLanded = world.testLanded(marsLander: nextLander)
				if isLanded {
					chromosome.ttl = Chromosome.maxTTL + nextLander.fuel
					chromosome.blackBox = nextLander
					return
				} else {
					let alternativeLander = world.simulate(marsLander: evolvingLander, action: Chromosome.neutralGene.action)
					if world.testLanded(marsLander: alternativeLander) {
						let prefixChromosome = chromosome.prefix(n: ttl-1)
						chromosome = prefixChromosome
						chromosome.ttl = Chromosome.maxTTL + alternativeLander.fuel
						chromosome.blackBox = alternativeLander
					} else {
						chromosome.ttl = ttl
						chromosome.blackBox = nextLander
					}
					return
				}
			}
		}
		log("Assertion acieed with chromosome: \(chromosome)")
		assert(false)
	}
}

extension Generation {

	func randomAction() -> Chromosome.Action {
		return Chromosome.Action(rotate: random[-90...90], power: random[0...4])
	}

	func generateMonoChromosome(with action: Chromosome.Action) -> Chromosome {
		return Chromosome(genes: [(action: action, duration: Chromosome.maxTTL)], ttl: -1, blackBox: nil)
	}

	mutating func populateToLimitWithRandom() {
		current.reserveCapacity(populationLimit * 2)
		for _ in current.count..<populationLimit {
			let newChromosome = Chromosome(
				genes: [(action: randomAction(), duration: Chromosome.maxTTL)],
				ttl: -1,
				blackBox: nil)
			current.append(newChromosome)
		}

		fitnessScored = false
	}

	mutating func addMutantsWithRandomTailOrHead() {
		let countBeforeMutation = current.count
		current.append(contentsOf: repeatElement(Chromosome(), count: 3*countBeforeMutation))
		for i in 0..<countBeforeMutation {
			let sample = current[i]

			current[countBeforeMutation+3*i+0] = makeMutant(sample, species: .HeadMutant)
			current[countBeforeMutation+3*i+1] = makeMutant(sample, species: .TailMutant)
			current[countBeforeMutation+3*i+2] = makeMutant(sample, species: .BodyMutant)
		}

		fitnessScored = false
	}

	enum MutantSpecies {
		case TailMutant
		case HeadMutant
		case BodyMutant
	}

	func makeMutant(_ sample: Chromosome, species: MutantSpecies) -> Chromosome {
		let mono = generateMonoChromosome(with: randomAction())

		guard case let ttl = UInt32(sample.ttl), sample.ttl > 2 else { return mono }
		let cut1 = random[1 ... Int(ttl/2)] //Int(1 + arc4random_uniform(ttl / 2))
		let cut2 = random[1 + Int(ttl/2) ... Int(ttl) - 1] //Int(1 + (ttl + 1) / 2 + arc4random_uniform((ttl - 2) / 2))

		switch( species ) {
		case .HeadMutant: return Chromosome.combine(head: mono, tail: sample, at: cut1)
		case .TailMutant: return Chromosome.combine(head: sample, tail: mono, at: cut2)
		case .BodyMutant: do {
			let withHead = Chromosome.combine(head: sample, tail: mono, at: cut1)
			let withHeadAndTail = Chromosome.combine(head: withHead, tail: sample, at: cut2)
			return withHeadAndTail
			}
		}
	}


	mutating func reducePopulation() {
		reducePopulation(to: populationLimit)
	}

	mutating func reducePopulation(to limit: Int) {
		guard current.count > limit else { return }

		evalTTL()
		sort()
		let newCurrent = sortedIndexes
			.prefix(limit)
			.map { current[$0] }
		current = newCurrent
		fitnessScored = false
	}

	mutating func bestChomosome() -> Chromosome {
		evalTTL()
		sort()
		return current[sortedIndexes[0]]
	}

	mutating func incrementAge(marsLander next: MarsLander) {
		lander = next
		for i in 0..<current.count {
			current[i].incrementAge()
		}

		fitnessScored = false
	}
}

extension Generation {

	mutating func evolution(cycles n: Int) {
		for _ in 0..<n {
			reducePopulation(to: populationLimit/5)
			addMutantsWithRandomTailOrHead()
			populateToLimitWithRandom()
		}
	}
}


let world = World.parseFromInput()
var action: MarsLander.Action!
var landerExpected: MarsLander!
var generation: Generation!
for turn in 0..<3000 {
	defer {
		print("\(action.rotate) \(action.power)")
		log("Expect Mars Lander: \(landerExpected!)")
	}

	var lander: MarsLander
	if feof(stdin) == 0 {
		lander = MarsLander(parseFromInput: readLine()!)!
		log("Read   Mars Lander: \(lander)")
	}
	else {
		lander = landerExpected
	}

	if let landerExpected = landerExpected {
		if lander - landerExpected < 2.0 {
			lander = landerExpected
		} else {
			log("Not close to expected, using data from input")
		}
	}

	let testAngle = { (angle: Int) in abs(lander.rotate - angle) < 90 ? true : false }

	// Apply DNA using fitness evaluation functions
	if turn == 0 {
		log("Calculating traectory \((lander.X, lander.Y)) -> \(world.target): ...")
		generation = Generation(world: world, lander: lander)

		let fitnessVX: Generation.ErrorFn = { (lander) in
			if (-20...20) ~= lander.hSpeed {
				return 0.0
			} else {
				return fabs(lander.hSpeed) - 20
			}
		}
		let fitnessVY: Generation.ErrorFn = { (lander) in
			if (-40...0) ~= lander.vSpeed {
				return 0.0
			} else {
				return fabs(lander.vSpeed) - 40
			}
		}
		let straightLanding: Generation.ErrorFn = { (lander) in
			if lander.rotate == 0 {
				return 0.0
			} else {
				return 1.0
			}
		}
		let radar: Generation.ErrorFn = { (lander) in
			let (x, y) = (lander.X, lander.Y)
			if (0...7000) ~= x && (0...3000) ~= y {
				return 0.0
			} else {
				return 1.0
			}
		}
		let targetX = world.target.x
		let deltaX: Generation.ErrorFn = { (lander: MarsLander) in
			let dist = fabs(targetX - lander.X)
			if dist < 500.0 {
				return 0.0
			} else {
				return Double(Int((dist - 500.0) / 10.0))
			}
		}
		let gluttony: Generation.ErrorFn = { (lander: MarsLander) in
			return Double(lander.fuel)
		}
		let deadman: Generation.ErrorFn = { (lander: MarsLander) in
			if lander.fuel <= 100 {
				return 0.0
			} else {
				return 1.0 / Double(lander.fuel - 100)
			}
		}

		generation.fitnessFunc.append((fn: fitnessVX, weight: 4.0))
		generation.fitnessFunc.append((fn: fitnessVY, weight: 4.0))
		generation.fitnessFunc.append((fn: straightLanding, weight: 0.5))
		generation.fitnessFunc.append((fn: radar, weight: 8.0))
		generation.fitnessFunc.append((fn: deltaX, weight: 10.0))
		generation.fitnessFunc.append((fn: gluttony, weight: 10.0))
		generation.fitnessFunc.append((fn: deadman, weight: 5.0))

		generation.populateToLimitWithRandom()
		generation.evolution(cycles: 10)
	}
	else {
		generation.incrementAge(marsLander: lander)
		generation.evalTTL()
		generation.evolution(cycles: 3)
		log("fitnessScore = \(generation.fitnessScore.prefix(upTo: 6).map{MarsLander.doubleFormatter($0)})")
	}

	let best = generation.bestChomosome()
	action = best.genes[0].action

	landerExpected = world.simulate(marsLander: lander, action: action)
}
